/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPartialVolumeModeller.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPartialVolumeModeller - convert an arbitrary dataset to a voxel representation
// with partial volume effects computed.
// .SECTION Description
// vtkPartialVolumeModeller is a filter that converts an arbitrary data set to a
// structured point (i.e., voxel) representation. It is very similar to
// vtkImplicitModeller, except that it doesn't record distance; instead it
// records the volume of the intersection of each voxel with the data set.
// .SECTION see also
// vtkImplicitModeller vtkVoxelModeller

#ifndef __vtkPartialVolumeModeller_h
#define __vtkPartialVolumeModeller_h

#include "vtkImageAlgorithm.h"

class vtkMultiThreader;
class vtkSimpleCriticalSection;

class VTK_ABI_EXPORT vtkPartialVolumeModeller : public vtkImageAlgorithm
{
public:
  vtkTypeMacro(vtkPartialVolumeModeller,vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Construct an instance of vtkPartialVolumeModeller with its sample dimensions
  // set to (50,50,50), and so that the model bounds are
  // automatically computed from its input. The maximum distance is set to
  // examine the whole grid. This could be made much faster, and probably
  // will be in the future.
  static vtkPartialVolumeModeller *New();

  // Description:
  // Compute the ModelBounds based on the input geometry.
  double ComputeModelBounds(double origin[3], double ar[3]);

  // Description:
  // Set the i-j-k dimensions on which to sample the distance function.
  // Default is (50, 50, 50)
  void SetSampleDimensions(int i, int j, int k);
  void SetSampleDimensions(int dim[3]);
  vtkGetVectorMacro(SampleDimensions,int,3);

  // Description:
  // Specify distance away from surface of input geometry to sample. Smaller
  // values make large increases in performance. Default is 1.0.
  vtkSetClampMacro(MaximumDistance,double,0.0,1.0);
  vtkGetMacro(MaximumDistance,double);

  // Description:
  // Specify the region in space to perform the voxelization.
  // Default is (0, 0, 0, 0, 0, 0)
  void SetModelBounds(double bounds[6]);
  void SetModelBounds(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
  vtkGetVectorMacro(ModelBounds,double,6);

    // Description:
  // Control the scalar type of the output image. The default is
  // VTK_DOUBLE.
  vtkSetMacro(OutputScalarType,int);
  void SetOutputScalarTypeToFloat(){this->SetOutputScalarType(VTK_FLOAT);};
  void SetOutputScalarTypeToDouble(){this->SetOutputScalarType(VTK_DOUBLE);};
  void SetOutputScalarTypeToInt(){this->SetOutputScalarType(VTK_INT);};
  void SetOutputScalarTypeToUnsignedInt()
    {this->SetOutputScalarType(VTK_UNSIGNED_INT);};
  void SetOutputScalarTypeToLong(){this->SetOutputScalarType(VTK_LONG);};
  void SetOutputScalarTypeToUnsignedLong()
    {this->SetOutputScalarType(VTK_UNSIGNED_LONG);};
  void SetOutputScalarTypeToShort(){this->SetOutputScalarType(VTK_SHORT);};
  void SetOutputScalarTypeToUnsignedShort()
    {this->SetOutputScalarType(VTK_UNSIGNED_SHORT);};
  void SetOutputScalarTypeToUnsignedChar()
    {this->SetOutputScalarType(VTK_UNSIGNED_CHAR);};
  void SetOutputScalarTypeToChar()
    {this->SetOutputScalarType(VTK_CHAR);};
  vtkGetMacro(OutputScalarType,int);

protected:
  vtkPartialVolumeModeller();
  ~vtkPartialVolumeModeller();

  virtual int RequestInformation(vtkInformation *,
                                 vtkInformationVector **,
                                 vtkInformationVector *);

    // see vtkAlgorithm for details
  virtual int RequestData(vtkInformation *request,
                          vtkInformationVector** inputVector,
                          vtkInformationVector* outputVector);

  // see algorithm for more info
  virtual int FillInputPortInformation(int port, vtkInformation* info);

  // Update a thread's progress. The parameter should range between
  // (0,1), and it represents the fraction of total work done by the
  // thread (not total work done by the filter).
  virtual void UpdateThreadProgress(double threadProgress);

  static VTK_THREAD_RETURN_TYPE ThreadedExecute( void *arg );

  vtkMultiThreader         *Threader;
  int                       NumberOfThreads;
  vtkSimpleCriticalSection *ProgressMutex;

  int    SampleDimensions[3];
  double MaximumDistance;
  double ModelBounds[6];
  int    OutputScalarType;

  // Keeps track of the total progress of the filter
  double TotalProgress;

private:
  vtkPartialVolumeModeller(const vtkPartialVolumeModeller&); // Not implemented
  void operator=(const vtkPartialVolumeModeller&); // Not implemented
};

#endif // __vtkPartialVolumeModeller_h
