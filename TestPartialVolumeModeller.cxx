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

#include <iostream>

#include <vtkCommand.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkImageShiftScale.h>
#include <vtkPNGWriter.h>
#include <vtkPoints.h>
#include <vtkRectilinearGridToTetrahedra.h>
#include <vtkTetra.h>
#include <vtkType.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVoxelModeller.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include "vtkPartialVolumeModeller.h"


class vtkProgressCommand : public vtkCommand
{
public:
  static vtkProgressCommand *New(){
    return new vtkProgressCommand;
  }

  virtual void Execute(vtkObject *caller, unsigned long, void *callData){
    double progress = *(static_cast<double*>(callData));
    fprintf(stderr, "\rFilter progress: %5.1f\n", 100.0 * progress);
    std::cerr.flush();
  }
};


int main(int argc, char* argv[])
{
  // Create a box partitioned into tetrahedra.
  double spacing = 5.0;
  double boxVoxelsCoveredX = 20;
  double boxVoxelsCoveredY = 30;
  double boxVoxelsCoveredZ = 15;

  double boxWidth = spacing * boxVoxelsCoveredX;
  double boxHeight = spacing * boxVoxelsCoveredY;
  double boxDepth = spacing * boxVoxelsCoveredZ;

  // Node numbering scheme
  //       7-----6
  //      /|    /|
  //     4-+---5 |
  //     | 3---+-2
  //     |/    |/
  //     0-----1
  double tetPoints[8][3] = {
    {0.0, 0.0, boxDepth},
    {boxWidth, 0.0, boxDepth},
    {boxWidth, 0.0, 0.0},
    {0.0, 0.0, 0.0},
    {0.0, boxHeight, boxDepth},
    {boxWidth, boxHeight, boxDepth},
    {boxWidth, boxHeight, 0.0},
    {0.0, boxHeight, 0.0}};

  vtkIdType tetIndices[6][4] = {
    {1, 3, 0, 4},
    {4, 1, 3, 5},
    {4, 7, 5, 3},
    {5, 7, 6, 3},
    {1, 3, 5, 2},
    {3, 2, 6, 5}
  };

  vtkUnstructuredGrid *grid = vtkUnstructuredGrid::New();

  vtkPoints* gridPoints = vtkPoints::New();

  for (int i = 0; i < 8; i++)
    {
    gridPoints->InsertNextPoint(tetPoints[i]);
    }
  grid->SetPoints(gridPoints);

  for (int i = 0; i < 6; i++)
    {
    grid->InsertNextCell(VTK_TETRA, 4, tetIndices[i]);
    }

  vtkXMLUnstructuredGridWriter *gridWriter = vtkXMLUnstructuredGridWriter::New();
  gridWriter->SetFileName("TetrahedralGrid.vtu");
  gridWriter->SetInput(grid);
  gridWriter->Write();

  // Set up partial volume modeller to create an image larger than the
  // extent of the box.
  double padding = 3 * spacing;
  double off = 0.0; // An offset applied to all modeller bounds
  double xmin = -padding + off;
  double xmax = boxWidth + padding + off;
  double ymin = -padding + off;
  double ymax = boxHeight + padding + off;
  double zmin = -padding + off;
  double zmax = boxDepth + padding + off;

  int imageDimsX = static_cast<int>((xmax - xmin) / spacing) + 1;
  int imageDimsY = static_cast<int>((ymax - ymin) / spacing) + 1;
  int imageDimsZ = static_cast<int>((zmax - zmin) / spacing) + 1;

  std::cout << "Creating image of dimensions (" << imageDimsX << ", " << imageDimsY << ", "
            << imageDimsZ << ")" << std::endl;

  vtkProgressCommand * pobserver = vtkProgressCommand::New();

  // Try out the vtkPartialVolumeModeller
  vtkPartialVolumeModeller *partialVolumeModeller = vtkPartialVolumeModeller::New();
  partialVolumeModeller->SetModelBounds(xmin, xmax, ymin, ymax, zmin, zmax);
  partialVolumeModeller->SetSampleDimensions(imageDimsX, imageDimsY, imageDimsZ);
  partialVolumeModeller->SetInput(grid);
  partialVolumeModeller->AddObserver(vtkCommand::ProgressEvent, pobserver);
  partialVolumeModeller->Update();
  pobserver->Delete();

  // Check out the results at the side edge of the box. Should be 0.25.
  double edgeValue = partialVolumeModeller->GetOutput()->
    GetScalarComponentAsDouble(3, 3, 4, 0);
  if (fabs(edgeValue - 0.25) > 1e-5)
    {
    std::cerr << "Expected value at voxel (3, 3, 4) should be 0.25, got "
              << edgeValue << std::endl;
    return EXIT_FAILURE;
    }

  // Check out the results at the corner of the box. Should be 0.125.
  double cornerValue = partialVolumeModeller->GetOutput()->
    GetScalarComponentAsDouble(3, 3, 3, 0);
  if (fabs(cornerValue - 0.125) > 1e-5)
    {
    std::cerr << "Expected value at voxel (3, 3, 3) to be 0.125, got "
              << cornerValue << std::endl;
    return EXIT_FAILURE;
    }

  vtkXMLImageDataWriter* imageWriter = vtkXMLImageDataWriter::New();
  imageWriter->SetFileName("PartialVolumeImageData.vti");
  imageWriter->SetInputConnection(partialVolumeModeller->GetOutputPort());
  imageWriter->Write();

  // Scale values to range necessary for image file formats.
  vtkImageShiftScale *shiftScale = vtkImageShiftScale::New();
  shiftScale->SetShift(0.0);
  shiftScale->SetScale(255.0);
  shiftScale->SetOutputScalarTypeToUnsignedChar();
  shiftScale->SetInputConnection(partialVolumeModeller->GetOutputPort());

  // Save as a PNG.
  vtkPNGWriter *writer = vtkPNGWriter::New();
  writer->SetFilePattern("PartialVolumeModellerImage%04d.png");
  writer->SetFileDimensionality(2);
  writer->SetInputConnection(shiftScale->GetOutputPort());
  writer->Write();

  // Try the vtkVoxelModeller for comparison
  vtkVoxelModeller *voxelModeller = vtkVoxelModeller::New();
  voxelModeller->SetModelBounds(xmin, xmax, ymin, ymax, zmin, zmax);
  voxelModeller->SetSampleDimensions(imageDimsX, imageDimsY, imageDimsZ);
  voxelModeller->SetScalarTypeToFloat();
  voxelModeller->SetInput(grid);
  voxelModeller->Update();

  imageWriter->SetFileName("VoxelModellerImageData.vti");
  imageWriter->SetInputConnection(partialVolumeModeller->GetOutputPort());
  imageWriter->Write();

  shiftScale->SetInputConnection(voxelModeller->GetOutputPort());

  writer->SetFilePattern("VoxelModellerImage%04d.png");
  writer->Write();

  writer->Delete();
  imageWriter->Delete();
  shiftScale->Delete();
  voxelModeller->Delete();

  partialVolumeModeller->Delete();

  return EXIT_SUCCESS;
}
