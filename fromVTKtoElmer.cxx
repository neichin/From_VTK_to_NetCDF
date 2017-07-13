#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <iostream>
#include "vtkXMLPUnstructuredGridReader.h"
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkPointSet.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkTriangle.h>
#include <vtkType.h>
#include "vtkConfigure.h"
#include "vtk_kwiml.h"
#include <vector>
#include <netcdfcpp.h>


using namespace std;

int isGrounded (vtkDoubleArray* data, vtkTriangle* triangle)
{
  double dataP1 = data->GetValue(triangle->GetPointIds()->GetId(0));
  double dataP2 = data->GetValue(triangle->GetPointIds()->GetId(1));
  double dataP3 = data->GetValue(triangle->GetPointIds()->GetId(2));

  if ((dataP1+dataP2+dataP3) == 3.0)
  {
	return 1;
  }
  else
  {
	return 0;
  }
}

double interpoleValue (vtkDoubleArray* data, vtkTriangle* triangle, double* weights)
{
  double interpoleVal;
  double dataP1 = data->GetValue(triangle->GetPointIds()->GetId(0));
  double dataP2 = data->GetValue(triangle->GetPointIds()->GetId(1));
  double dataP3 = data->GetValue(triangle->GetPointIds()->GetId(2));

  interpoleVal = dataP1*weights[0] + dataP2*weights[1] + dataP3*weights[2];


  return interpoleVal;
}

vtkSmartPointer<vtkTriangle> findTriangle (vtkUnstructuredGrid *output, double  point[3], double *pcoords, double *weights, double tol2=5000.0, int subid=0)
{
  vtkSmartPointer<vtkTriangle> tri = vtkSmartPointer<vtkTriangle>::New();
  vtkSmartPointer<vtkCell> cell = output->FindAndGetCell(point,tri,5,tol2,subid,pcoords,weights);
  vtkSmartPointer<vtkTriangle> triangle = vtkTriangle::SafeDownCast(cell);
  
  return triangle;
}

static const int NC_ERR = 2;
static const NcToken xDimNAME = "x";
static const NcToken yDimNAME = "y";

static const NcToken xVarNAME = "x";
static const NcToken yVarNAME = "y";
static const NcToken draftVarNAME = "isf_draft";

int main(int argc, char* argv[])
{
  string filename;
  string filenameNC;
  string filenameOutput;
  string filenameBathy;
  double xmin = 350000.0;
  double xmax = 640000.0;
  double ymin = 0.0;
  double ymax = 80000.0;
  //////////// ARGS PARSING
  if (argc < 5)
  {
	cout << "Error: at least 5 arguments are required " << endl;
	return 0;
  }
  else 
  {
	filename = argv[1];
 	filenameNC = argv[2];
	filenameOutput = argv[4];
	filenameBathy = argv[3];
  }

  if (argc > 5 && argc != 9)
  {
	cout << "Error: 4 arguments required for domain bounds " << endl;
        return 0;
  }
  else if (argc > 5 )
  {
	xmin = atof(argv[5]);
	xmax = atof(argv[6]);
	ymin = atof(argv[7]);
	ymax = atof(argv[8]);
  }

  cout << "Interpolating Elmer to NEMO grid" << endl;
  cout << "           " << endl;
  cout << "Elmer grid is found in; " << endl;
  cout << "        " << filename << endl;
 
  cout << xmin <<" "<< xmax << " " << ymin << " " << ymax << " " << endl;
  /////////////////////////////////////////// 
  //OPEN AND READ VTK FILE
  vtkSmartPointer<vtkXMLPUnstructuredGridReader> source = vtkSmartPointer<vtkXMLPUnstructuredGridReader>::New();
  source->SetFileName(filename.c_str());
  source->Update();
  vtkSmartPointer<vtkUnstructuredGrid> output = source->GetOutput();
  vtkSmartPointer<vtkPointData> pointsData = output->GetPointData();

  vtkSmartPointer<vtkDataArray> arrayDraft = pointsData->GetScalars("zb");
  vtkSmartPointer<vtkDoubleArray> arrayDoubleDraft = vtkDoubleArray::SafeDownCast(arrayDraft);


  vtkSmartPointer<vtkDataArray> arrayGM = pointsData->GetScalars("groundedmask");
  vtkSmartPointer<vtkDoubleArray> arrayDoubleGM = vtkDoubleArray::SafeDownCast(arrayGM);
  ////////////////////////
  //////////////////////////////////////////
  
  cout << "OPEN NETCDF....." << endl;

  /////////////////////////////////////////
  //OPEN AND READ NETCDF GRID
  NcError err(NcError::silent_nonfatal);
  NcFile dataFile(filenameNC.c_str(), NcFile::ReadOnly);

  NcFile dataFileBathy(filenameBathy.c_str(), NcFile::ReadOnly);
  
  NcVar *bathy, *x, *y, *draft; 
  //GET DIMENSIONS
  NcDim *xdim = dataFile.get_dim(0);
  NcDim *ydim = dataFile.get_dim(1);
  long nx = xdim->size();
  long ny = ydim->size();

  cout << "NetCDF dims " << nx << " x " << ny << endl;

  // READ X DATA
  x = dataFile.get_var(xVarNAME);
  float xData[nx];
  x->get(xData,nx);

  // READ Y DATA
  y = dataFile.get_var(yVarNAME);
  float yData[nx];
  y->get(yData,ny);

  // READ BATHY DATA
  bathy = dataFileBathy.get_var("Bathymetry_isf");
  float BathyData[ny][nx];
  bathy->get(&BathyData[0][0],ny,nx);

  /////////////////////////////////////////////
  //CREATE ARRAY FOR NEW DRAFT
  float newDraft[ny][nx];
  vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
  double point[3];
  double pcoords[3];
  double weights[3];

  for(int i=0; i<nx ; i++)
  {
    for(int j=0; j<ny ; j++)
    {
	point[0]=xData[i];
	point[1]=yData[j];
	point[3]=0.0;
	if ( point[1]< ymin || point[1]> ymax )
	{
	  newDraft[j][i] = 0.0 ;
	  continue;
	}
        if ( point[0]> xmax )
        //if ( point[0]< xmin || point[0]> xmax )
        {
	  newDraft[j][i] = 0.0;
          continue;
        }
	triangle = findTriangle(output, point, pcoords, weights);
	if (isGrounded(arrayDoubleGM, triangle))
	{
	  newDraft[j][i] = BathyData[j][i];
	}else{
	  newDraft[j][i] = -1.0 * interpoleValue(arrayDoubleDraft, triangle, weights);
	}
    }
  }
  ////////////////////////////////////////////      
  ////
  //NcFile dataFileNew(filenameNC.c_str(), NcFile::Replace);
  //
  cout << "Reading Netcdf" << endl;

  NcFile dataFileNew(filenameOutput.c_str(), NcFile::Replace);
  NcDim *xDim = dataFileNew.add_dim("x", nx);
  NcDim *yDim = dataFileNew.add_dim("y", ny);
  NcVar *xNew = dataFileNew.add_var("x", ncFloat, xDim);
  NcVar *yNew = dataFileNew.add_var("y", ncFloat, yDim);
  NcVar *draftNew = dataFileNew.add_var("isf_draft", ncFloat, yDim, xDim);
  xNew->put(xData,nx);
  yNew->put(yData,ny);
  draftNew->put(&newDraft[0][0],ny,nx);








  return EXIT_SUCCESS;
}
