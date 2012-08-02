/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPartialVolumeModeller.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkPartialVolumeModeller.h"

#include "vtkBoxClipDataSet.h"
#include "vtkCell.h"
#include "vtkCellLocator.h"
#include "vtkCriticalSection.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkGenericCell.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiThreader.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStructuredGrid.h"
#include "vtkTetra.h"
#include "vtkUnstructuredGrid.h"

#include <math.h>

vtkStandardNewMacro(vtkPartialVolumeModeller);

struct vtkPartialVolumeModellerThreadInfo
{
  vtkPartialVolumeModeller *Modeller;
  vtkDataSet               **Input;
};

// Construct an instance of vtkPartialVolumeModeller with its sample dimensions
// set to (50,50,50), and so that the model bounds are
// automatically computed from its input. The maximum distance is set to
// examine the whole grid. This could be made much faster, and probably
// will be in the future.
vtkPartialVolumeModeller::vtkPartialVolumeModeller()
{
  this->MaximumDistance = 1.0;

  this->ModelBounds[0] = 0.0;
  this->ModelBounds[1] = 0.0;
  this->ModelBounds[2] = 0.0;
  this->ModelBounds[3] = 0.0;
  this->ModelBounds[4] = 0.0;
  this->ModelBounds[5] = 0.0;

  this->SampleDimensions[0] = 50;
  this->SampleDimensions[1] = 50;
  this->SampleDimensions[2] = 50;

  this->OutputScalarType = VTK_DOUBLE;

  this->Threader        = vtkMultiThreader::New();
  this->NumberOfThreads = this->Threader->GetNumberOfThreads();
  this->ProgressMutex = vtkSimpleCriticalSection::New();
}

//----------------------------------------------------------------------------
vtkPartialVolumeModeller::~vtkPartialVolumeModeller()
{
  if (this->Threader)
    {
    this->Threader->Delete();
    }

  if (this->ProgressMutex)
    {
    this->ProgressMutex->Delete();
    }
}

//----------------------------------------------------------------------------
// Specify the position in space to perform the voxelization.
void vtkPartialVolumeModeller::SetModelBounds(double bounds[6])
{
  this->SetModelBounds(bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);
}

//----------------------------------------------------------------------------
void vtkPartialVolumeModeller::SetModelBounds(double xmin, double xmax, double ymin,
                                      double ymax, double zmin, double zmax)
{
  if (this->ModelBounds[0] != xmin || this->ModelBounds[1] != xmax ||
      this->ModelBounds[2] != ymin || this->ModelBounds[3] != ymax ||
      this->ModelBounds[4] != zmin || this->ModelBounds[5] != zmax )
    {
    this->Modified();
    this->ModelBounds[0] = xmin;
    this->ModelBounds[1] = xmax;
    this->ModelBounds[2] = ymin;
    this->ModelBounds[3] = ymax;
    this->ModelBounds[4] = zmin;
    this->ModelBounds[5] = zmax;
    }
}

//----------------------------------------------------------------------------
int vtkPartialVolumeModeller::RequestInformation (
  vtkInformation * vtkNotUsed(request),
  vtkInformationVector ** vtkNotUsed( inputVector ),
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  int i;
  double ar[3], origin[3];

  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
               0, this->SampleDimensions[0]-1,
               0, this->SampleDimensions[1]-1,
               0, this->SampleDimensions[2]-1);

  for (i=0; i < 3; i++)
    {
    origin[i] = this->ModelBounds[2*i];
    if ( this->SampleDimensions[i] <= 1 )
      {
      ar[i] = 1;
      }
    else
      {
      ar[i] = (this->ModelBounds[2*i+1] - this->ModelBounds[2*i])
              / (this->SampleDimensions[i] - 1);
      }
    }
  outInfo->Set(vtkDataObject::ORIGIN(),origin,3);
  outInfo->Set(vtkDataObject::SPACING(),ar,3);

  vtkDataObject::SetPointDataActiveScalarInfo(outInfo, this->OutputScalarType, 1);
  return 1;
}

//----------------------------------------------------------------------------
VTK_THREAD_RETURN_TYPE vtkPartialVolumeModeller::ThreadedExecute( void *arg )
{
  int threadId = ((vtkMultiThreader::ThreadInfo *)(arg))->ThreadID;
  int threadCount = ((vtkMultiThreader::ThreadInfo *)(arg))->NumberOfThreads;
  vtkPartialVolumeModellerThreadInfo *userData = (vtkPartialVolumeModellerThreadInfo *)
    (((vtkMultiThreader::ThreadInfo *)(arg))->UserData);

  vtkDataSet *input = userData->Input[threadId];

  // Extract the grid boundaries
  vtkDataSetSurfaceFilter *surfaceFilter = vtkDataSetSurfaceFilter::New();
  surfaceFilter->SetInput(input);
  surfaceFilter->Update();

  // Set up a locator
  vtkCellLocator *locator = vtkCellLocator::New();
  locator->SetDataSet(surfaceFilter->GetOutput());
  locator->AutomaticOn();
  locator->SetNumberOfCellsPerBucket(1);
  locator->CacheCellBoundsOn();
  locator->BuildLocator();

  // Generic cell to use when locating cells
  vtkGenericCell *cell = vtkGenericCell::New();

  vtkImageData *output = userData->Modeller->GetOutput();
  double *spacing = output->GetSpacing();
  double *origin = output->GetOrigin();

  int *sampleDimensions = userData->Modeller->GetSampleDimensions();
  if (!output->GetPointData()->GetScalars())
    {
    vtkGenericWarningMacro("No output scalars defined.");
    surfaceFilter->Delete();
    locator->Delete();
    cell->Delete();

    return VTK_THREAD_RETURN_VALUE;
    }

  // break up into slabs based on threadId and threadCount
  int slabSize = sampleDimensions[2] / threadCount;
  int slabSizeExtra = slabSize + 1; // For slabs with an extra layer
  int numThreadsWithExtraLayer = sampleDimensions[2] % threadCount;

  int slabMin, slabMax;
  if (threadId < numThreadsWithExtraLayer)
    {
    slabMin = threadId * slabSizeExtra;
    slabMax = slabMin + slabSizeExtra - 1;
    }
  else
    {
    slabMin = numThreadsWithExtraLayer * slabSizeExtra +
      (threadId - numThreadsWithExtraLayer) * slabSize;
    slabMax = slabMin + slabSize - 1;
    }
  if (slabMin >= sampleDimensions[2])
    {
    return VTK_THREAD_RETURN_VALUE;
    }

  double *weights = new double[input->GetMaxCellSize()];
  vtkDataArray *newScalars = output->GetPointData()->GetScalars();

  //
  // Voxel widths are 1/2 the height, width, length of a voxel
  //
  double voxelHalfWidth[3];
  for (int i=0; i < 3; i++)
    {
    voxelHalfWidth[i] = spacing[i] / 2.0;
    }

  double voxelRadius = sqrt(voxelHalfWidth[0]*voxelHalfWidth[0] +
                            voxelHalfWidth[1]*voxelHalfWidth[1] +
                            voxelHalfWidth[2]*voxelHalfWidth[2]);

  // Set up the box clipping filter
  vtkBoxClipDataSet *clipper = vtkBoxClipDataSet::New();
  clipper->SetInput(input);
  clipper->SetOrientation(0);

  // Compute the volume of a filled voxel.
  double fullVoxelVolume = spacing[0]*spacing[1]*spacing[2];

  //
  // Traverse all voxels, computing partial volume intersecton on all points.
  //
  int jkFactor = sampleDimensions[0]*sampleDimensions[1];
  double threadTotalVoxels = static_cast< double >( (slabMax - slabMin + 1) * jkFactor );
  int numThreads = userData->Modeller->Threader->GetNumberOfThreads();
  double voxelProgressWeight = 1.0 / threadTotalVoxels;
  double voxelPoint[3];
  int count = 0;
  for (int k = slabMin; k <= slabMax; k++)
    {
    voxelPoint[2] = static_cast<double>(k)*spacing[2] + origin[2];
    double zmin = voxelPoint[2] - voxelHalfWidth[2];
    double zmax = voxelPoint[2] + voxelHalfWidth[2];
    for (int j = 0; j < sampleDimensions[1]; j++)
      {
      voxelPoint[1] = static_cast<double>(j)*spacing[1] + origin[1];
      double ymin = voxelPoint[1] - voxelHalfWidth[1];
      double ymax = voxelPoint[1] + voxelHalfWidth[1];
      for (int i = 0; i < sampleDimensions[0]; i++)
        {
        voxelPoint[0] = static_cast<double>(i)*spacing[0] + origin[0];
        double xmin = voxelPoint[0] - voxelHalfWidth[0];
        double xmax = voxelPoint[0] + voxelHalfWidth[0];

        // Look for the closest point on the surface. If it is outside the voxel,
        // we can skip the expensive box clipping.
        double closestPoint[3];
        vtkIdType cellId;
        int subId;
        double dist2;
        int cellFound =
          locator->FindClosestPointWithinRadius(voxelPoint, voxelRadius,
                                                closestPoint, cell, cellId,
                                                subId, dist2);

        bool boundaryIntersectsVoxel = cellFound == 1 &&
          closestPoint[0] >= xmin && closestPoint[0] <= xmax &&
          closestPoint[1] >= ymin && closestPoint[1] <= ymax &&
          closestPoint[2] >= zmin && closestPoint[2] <= zmax;

        double volume = 0.0;
        if (boundaryIntersectsVoxel)
          {
          // We need to clip with the voxel bounds and compute the volume of the intersection
          // Update the box dimensions in the clipper
          clipper->SetBoxClip(voxelPoint[0] - voxelHalfWidth[0], voxelPoint[0] + voxelHalfWidth[0],
                              voxelPoint[1] - voxelHalfWidth[1], voxelPoint[1] + voxelHalfWidth[1],
                              voxelPoint[2] - voxelHalfWidth[2], voxelPoint[2] + voxelHalfWidth[2]);
          clipper->Update();

          // Get the output tetrahedra from the clipper and compute the sum of their volumes
          vtkUnstructuredGrid* intersection = clipper->GetOutput();
          int numCells = intersection->GetNumberOfCells();
          for (int cellNum = 0; cellNum < numCells; cellNum++)
            {
            vtkCell *cell = intersection->GetCell(cellNum);
            if (cell->GetCellType() == VTK_TETRA)
              {
              vtkTetra* tet = vtkTetra::SafeDownCast(cell);
              vtkPoints* pts = tet->GetPoints();
              double p0[3], p1[3], p2[3], p3[3];
              pts->GetPoint(0, p0);
              pts->GetPoint(1, p1);
              pts->GetPoint(2, p2);
              pts->GetPoint(3, p3);
              double cellVolume = vtkTetra::ComputeVolume(p0, p1, p2, p3);

              // Should really check for negative volume elements here, but it's not
              // clear what the right course of action is in that case.
              volume += cellVolume;
              }
            }
          }
        else // no boundary intersection with voxel
          {
          // The voxel is either completely inside the grid or completely outside the grid.
          // We need to determine which.
          double pcoords[3];
          cellId = input->FindCell(voxelPoint, NULL, 0, 1e-5, subId, pcoords, weights);
          if (cellId == -1)
            {
            volume = 0.0;
            }
          else
            {
            volume = fullVoxelVolume;
            }
          }

        int idx = jkFactor*k + sampleDimensions[0]*j + i;
        newScalars->SetComponent(idx, 0, volume / fullVoxelVolume);

        if (count == 50)
          {
          userData->Modeller->UpdateThreadProgress(voxelProgressWeight*count);
          count = 0;

          if (threadId == 0)
            {
            userData->Modeller->UpdateProgress( userData->Modeller->TotalProgress );
            }
          }
        ++count;
        }
      }
    }

  // Report the remnants
  userData->Modeller->UpdateThreadProgress(voxelProgressWeight*count);
  if (threadId == 0)
    {
    userData->Modeller->UpdateProgress( userData->Modeller->TotalProgress );
    }

  clipper->Delete();
  surfaceFilter->Delete();
  locator->Delete();
  cell->Delete();

  delete [] weights;

  return VTK_THREAD_RETURN_VALUE;
}

//----------------------------------------------------------------------------
int vtkPartialVolumeModeller::RequestData(
  vtkInformation* vtkNotUsed( request ),
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // get the input
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkDataSet *input = vtkDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Check that the input type is supported
  int dataObjectType = input->GetDataObjectType();
  if (dataObjectType != VTK_STRUCTURED_GRID &&
      dataObjectType != VTK_UNSTRUCTURED_GRID &&
      dataObjectType != VTK_RECTILINEAR_GRID)
    {
    vtkErrorMacro(<< "vtkPartialVolumeModeller expects an input data set of type "
                  << "vtkStructuredGrid, vtkUnstructuredGrid, or vtkRectilinearGrid.");
    return 0;
    }

  // get the output
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkImageData *output = vtkImageData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // We need to allocate our own scalars since we are overriding
  // the superclasses "Execute()" method.
  output->SetExtent(output->GetWholeExtent());
  output->AllocateScalars();

  double origin[3], spacing[3];
  double maxDistance = this->ComputeModelBounds(origin, spacing);
  outInfo->Set(vtkDataObject::ORIGIN(), origin, 3);

  vtkPartialVolumeModellerThreadInfo info;
  info.Modeller = this;

  // Deep copy the data set to avoid a race condition. We could also
  // split the mesh into slabs to reduce the amount of work each thread
  // has to do, but the splitting needs to be done in one thread,
  // so doing so is not a clear win.
  info.Input = new vtkDataSet*[this->Threader->GetNumberOfThreads()];
  for (int threadId = 0; threadId < this->Threader->GetNumberOfThreads(); threadId++)
    {
    switch (dataObjectType)
      {
      case VTK_STRUCTURED_GRID:
        info.Input[threadId] = vtkStructuredGrid::New();
        break;

      case VTK_UNSTRUCTURED_GRID:
        info.Input[threadId] = vtkUnstructuredGrid::New();
        break;

      case VTK_RECTILINEAR_GRID:
        info.Input[threadId] = vtkRectilinearGrid::New();
        break;
      }
    info.Input[threadId]->DeepCopy(input);
    }

  // Set the number of threads to use,
  // then set the execution method and do it.
  // If the number of threads is greater than the image dimension
  // along the splitting axis (z), reduce the number of threads so
  // that each one gets a single-layer slab
  if ( this->NumberOfThreads > this->SampleDimensions[2] )
    {
    this->Threader->SetNumberOfThreads( this->SampleDimensions[2] );
    }
  else
    {
    this->Threader->SetNumberOfThreads( this->NumberOfThreads );
    }
  this->Threader->SetSingleMethod( vtkPartialVolumeModeller::ThreadedExecute,
    (void *)&info);
  this->TotalProgress = 0.0;
  this->Threader->SingleMethodExecute();

  // Clean up.
  for (int threadId = 0; threadId < this->Threader->GetNumberOfThreads(); threadId++)
    {
    info.Input[threadId]->Delete();
    }
  delete[] info.Input;

  return 1;
}

//----------------------------------------------------------------------------
// Compute the ModelBounds based on the input geometry.
double vtkPartialVolumeModeller::ComputeModelBounds(double origin[3],
                                                    double spacing[3])
{
  double *bounds, maxDist;
  int i, adjustBounds=0;

  // compute model bounds if not set previously
  if ( this->ModelBounds[0] >= this->ModelBounds[1] ||
       this->ModelBounds[2] >= this->ModelBounds[3] ||
       this->ModelBounds[4] >= this->ModelBounds[5] )
    {
    adjustBounds = 1;
    vtkDataSet *ds = vtkDataSet::SafeDownCast(this->GetInput());
    // ds better be non null otherwise something is very wrong here
    bounds = ds->GetBounds();
    }
  else
    {
    bounds = this->ModelBounds;
    }

  for (maxDist=0.0, i=0; i<3; i++)
    {
    if ( (bounds[2*i+1] - bounds[2*i]) > maxDist )
      {
      maxDist = bounds[2*i+1] - bounds[2*i];
      }
    }
  maxDist *= this->MaximumDistance;

  // adjust bounds so model fits strictly inside (only if not set previously)
  if ( adjustBounds )
    {
    for (i=0; i<3; i++)
      {
      this->ModelBounds[2*i] = bounds[2*i] - maxDist;
      this->ModelBounds[2*i+1] = bounds[2*i+1] + maxDist;
      }
    }

  // Set volume origin and data spacing
  for (i=0; i<3; i++)
    {
    origin[i] = this->ModelBounds[2*i];
    spacing[i] = (this->ModelBounds[2*i+1] - this->ModelBounds[2*i])/
      (this->SampleDimensions[i] - 1);
    }

  return maxDist;
}

//----------------------------------------------------------------------------
// Set the i-j-k dimensions on which to sample the distance function.
void vtkPartialVolumeModeller::SetSampleDimensions(int i, int j, int k)
{
  int dim[3];

  dim[0] = i;
  dim[1] = j;
  dim[2] = k;

  this->SetSampleDimensions(dim);
}

//----------------------------------------------------------------------------
void vtkPartialVolumeModeller::SetSampleDimensions(int dim[3])
{
  int dataDim, i;

  vtkDebugMacro(<< " setting SampleDimensions to (" << dim[0] << "," << dim[1]
                << "," << dim[2] << ")");

  if ( dim[0] != this->SampleDimensions[0] ||
       dim[1] != this->SampleDimensions[1] ||
       dim[2] != this->SampleDimensions[2] )
    {
    if ( dim[0]<1 || dim[1]<1 || dim[2]<1 )
      {
      vtkErrorMacro (<< "Bad Sample Dimensions, retaining previous values");
      return;
      }
    for (dataDim=0, i=0; i<3 ; i++)
      {
      if (dim[i] > 1)
        {
        dataDim++;
        }
      }
    if ( dataDim  < 3 )
      {
      vtkErrorMacro(<<"Sample dimensions must define a volume!");
      return;
      }

    for ( i=0; i<3; i++)
      {
      this->SampleDimensions[i] = dim[i];
      }
    this->Modified();
    }
}

//----------------------------------------------------------------------------
int vtkPartialVolumeModeller::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}

//----------------------------------------------------------------------------
void vtkPartialVolumeModeller::UpdateThreadProgress(double threadProgress)
{
  this->ProgressMutex->Lock();

  // This doesn't exactly represent the fraction of work contributed
  // by the thread to the total problem, but it's good enough.
  this->TotalProgress += threadProgress / static_cast<double>(this->Threader->GetNumberOfThreads());

  this->ProgressMutex->Unlock();
}

//----------------------------------------------------------------------------
void vtkPartialVolumeModeller::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Maximum Distance: " << this->MaximumDistance << "\n";
  os << indent << "Sample Dimensions: (" << this->SampleDimensions[0] << ", "
               << this->SampleDimensions[1] << ", "
               << this->SampleDimensions[2] << ")\n";
  os << indent << "Model Bounds: \n";
  os << indent << "  Xmin,Xmax: (" << this->ModelBounds[0] << ", "
     << this->ModelBounds[1] << ")\n";
  os << indent << "  Ymin,Ymax: (" << this->ModelBounds[2] << ", "
     << this->ModelBounds[3] << ")\n";
  os << indent << "  Zmin,Zmax: (" << this->ModelBounds[4] << ", "
     << this->ModelBounds[5] << ")\n";
  os << indent << "OutputScalarType: " << this->OutputScalarType << endl;
}
