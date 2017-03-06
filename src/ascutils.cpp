/*
 *  ascutils.cpp
 *  ESRI ASCII grid utility functions
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/GRID
 *
 */

#include "ascutils.h"
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <string.h>
#include <math.h>
#include <fstream>
#include <streambuf>
#include <sstream>

//---------------------------------------------------------------------------------------------------------------

void asc_grid::expandstorages(asc_header hdr)
{
	int nrows=hdr.nrows;
	int ncols=hdr.ncols;
	sdata = new double [nrows*ncols];
	data = new double * [nrows];
	for(int r=0; r<nrows; r++){
		data[r] = sdata + r * ncols;
	}
}

//---------------------------------------------------------------------------------------------------------------

asc_grid::asc_grid (std::string filename)
{
	loadFromFile(filename);
}

//---------------------------------------------------------------------------------------------------------------

asc_grid::asc_grid (asc_grid * sample, bool takedata)
{
	hdr=sample->header();
	expandstorages(hdr);
	if(takedata){
		copyData(sample->data);
	}
	if(takedata){
		integerdata=sample->integerdata;
	}
}

//---------------------------------------------------------------------------------------------------------------

asc_grid::asc_grid (asc_header hdr1, double ** data1)
{
	hdr=hdr1;
	expandstorages(hdr);
	copyData(data1);
	integerdata=false;
}

//---------------------------------------------------------------------------------------------------------------

asc_grid::asc_grid ()
{
	//empty grid
	hdr.nrows = 0;
	hdr.ncols = 0;
	hdr.xll = 0;
	hdr.yll = 0;
	hdr.cellsize = 0;
	hdr.nodata_value = 0;
	
	sdata = NULL;
	data = NULL;
	integerdata=false;
}

//---------------------------------------------------------------------------------------------------------------

asc_grid::asc_grid (asc_header hdr1, double * sdata1)
{
	hdr=hdr1;
	expandstorages(hdr);
	if(sdata1){
		memcpy(sdata,sdata1,hdr.nrows*hdr.ncols*sizeof(double));
	}
}

//---------------------------------------------------------------------------------------------------------------

void asc_grid::copyData(double ** data1)
{
	if(data && data1){
		for(int r=0; r<hdr.nrows; r++){
			for(int c=0; c<hdr.ncols; c++){
				data[r][c]=data1[r][c];
			}
		}
	}
}

//---------------------------------------------------------------------------------------------------------------

double asc_grid::operator[](unsigned int index)
{
	double result=0.0;
	if(sdata && index<hdr.nrows*hdr.ncols){
		result=sdata[index];
	}
	return result;
}

//---------------------------------------------------------------------------------------------------------------

asc_grid& asc_grid::operator= (asc_grid const& g)
{
	if (this == &g){
		return *this;   // Gracefully handle self assignment
	}

	//delocate memory for old data
	if(sdata){
		delete [] sdata;
	}
	if(data){
		delete [] data;
	}
	
	//copy header 
	hdr = g.hdr;
	integerdata=g.integerdata;

	id=g.id;
	kind=g.kind;
	displayname=g.displayname;
	
	expandstorages(hdr);
	copyData(g->data)
	
	return *this;
} 

//---------------------------------------------------------------------------------------------------------------

void asc_grid::loadFromFile(std::string filename)
{
	integerdata=true;
	double foo;
	FILE * f=openasc_r(filename, hdr);
	if(f){
		expandstorages(hdr);
		if(data){
			for(int r=0; r<hdr.nrows; r++){
				for(int c=0; c<hdr.ncols; c++){
					double val;
					if(fscanf(f,"%lf",&val)!=1){
						printf("[Error]: Data error in %s at row #%d column #%d.\n",filename.c_str(),r,c);
						return;
					}
					else{
						if(val!=hdr.nodata_value && modf(val, &foo)!=0.0){
							integerdata=false;
						}
						data[r][c]=val;
					}
				}
			}
		}
		fclose(f);
	}
	else{
		printf("[Error]: failed to open %s for reading.\n",filename.c_str());
		return;
	}
}

//---------------------------------------------------------------------------------------------------------------

double asc_grid::min()
{
	double result=0.0;
	bool firstval=true;
	if(data){
		for(int r=0; r<hdr.nrows; r++){
			for(int c=0; c<hdr.ncols; c++){
				double val=data[r][c];
				if(val!=hdr.nodata_value){
					if(firstval){
						result=val;
						firstval=false;
					}
					else{
						if(val<result){
							result=val;
						}
					}
				}
			}
		}
	}
	return result;
}

//---------------------------------------------------------------------------------------------------------------

void asc_grid::setValue(double val)
{
	if(data){
		for(int r=0; r<hdr.nrows; r++){
			for(int c=0; c<hdr.ncols; c++){
				data[r][c]=val;
			}
		}
	}
}

//---------------------------------------------------------------------------------------------------------------

void asc_grid::fillClassified(std::map<int, double> values, asc_grid * class_map)
{
	//header check first
	if(!class_map || hdr.nrows!=class_map->hdr.nrows || hdr.ncols!=class_map->hdr.ncols || values.size()==0){
		return;
	}	
	
	//do the reclassification (acceptable for a few classes, uses simple caching)
	int class_nodata=class_map->hdr.nodata_value;
	double my_nodata=hdr.nodata_value; 
	int prev_class_val=class_nodata;
	double prev_value=my_nodata;
	double act_value;
	
	if(data && class_map->data){
		for(int r=0; r<hdr.nrows; r++){
			for(int c=0; c<hdr.ncols; c++){
				int class_val = class_map->data[r][c];
				if(class_val!=class_nodata){
					if(class_val==prev_class_val){
						act_value=prev_value;
					}
					else{
						act_value=values[class_val];
					}
				}
				else{
					act_value=my_nodata;
				}
				prev_class_val=class_val;
				prev_value=act_value;
				data[r][c]=act_value;
			}
		}
	}
}

//---------------------------------------------------------------------------------------------------------------

void asc_grid::fillClassified(std::map<int, double *> values, asc_grid * class_map)
{
	//header check first
	if(!class_map || hdr.nrows!=class_map->hdr.nrows || hdr.ncols!=class_map->hdr.ncols || values.size()==0){
		return;
	}	
	
	//do the reclassification (acceptable for a few classes, uses simple caching)
	int class_nodata=class_map->hdr.nodata_value;
	double my_nodata=hdr.nodata_value; 
	int prev_class_val=class_nodata;
	double prev_value=my_nodata;
	double act_value;
	
	if(data && class_map->data){
		for(int r=0; r<hdr.nrows; r++){
			for(int c=0; c<hdr.ncols; c++){
				int class_val = class_map->data[r][c];
				if(class_val!=class_nodata){
					if(class_val==prev_class_val){
						act_value=prev_value;
					}
					else{
						act_value=*(values[class_val]);	//here is the only difference with *()
					}
				}
				else{
					act_value=my_nodata;
				}
				prev_class_val=class_val;
				prev_value=act_value;
				data[r][c]=act_value;
			}
		}
	}
}

//---------------------------------------------------------------------------------------------------------------

double * asc_grid::ptrToCell(int r, int c)
{
	if(r>=0 && r<hdr.nrows && c>=0 && c<hdr.ncols){
		return &(data[r][c]);
	}
	return NULL;
}

//---------------------------------------------------------------------------------------------------------------

double * asc_grid::ptrToCoord(double x, double y)
{	
	double xul = hdr.xll;
	double yul = hdr.yll + hdr.nrows * hdr.cellsize;
	int r= hdr.nrows - (yul - y) / hdr.cellsize;
	int c= (x - xul) / hdr.cellsize;
	return ptrToCell(r,c);
}
	
//---------------------------------------------------------------------------------------------------------------

double asc_grid::max()
{
	double result=0.0;
	bool firstval=true;
	if(data){
		for(int r=0; r<hdr.nrows; r++){
			for(int c=0; c<hdr.ncols; c++){
				double val=data[r][c];
				if(val!=hdr.nodata_value){
					if(firstval){
						result=val;
						firstval=false;
					}
					else{
						if(val>result){
							result=val;
						}
					}
				}
			}
		}
	}
	return result;
}

//---------------------------------------------------------------------------------------------------------------

void asc_grid::addValuesFrom(asc_grid * grd)
{
	if(!grd || hdr.nrows!=grd->hdr.nrows || hdr.ncols!=grd->hdr.ncols){
		return;
	}
	
	double ** sdata2 = grd->sdata;
	int numdata = hdr.nrows * hdr.ncols;
	for(int w=0; w<numdata; w++){
		if(sdata[w]!=hdr.nodata_value && sdata2[w]!=grd->hdr.nodata_value){
			sdata[w]+=sdata2[w];
		}
		else{
			sdata[w]=hdr.nodata_value;
		}
	}
}
	
//---------------------------------------------------------------------------------------------------------------

asc_grid::~asc_grid()
{
	if(sdata){
		delete [] sdata;
		//collapseGrid(data, hdr);
	}
	if(data){
		delete [] data;
	}
}

//---------------------------------------------------------------------------------------------------------------

void asc_grid::saveToFile(std::string filename)
{
	FILE * f=openasc_w(filename, hdr);
	if(f){
		if(data){
			for(int r=0; r<hdr.nrows; r++){
				for(int c=0; c<hdr.ncols; c++){
					fprintf(f,"%g ",data[r][c]);
				}
			}
			fprintf(f,"\n");
		}
		fclose(f);
	}
	else{
		printf("[Error]: failed to open %s for writing.\n",filename.c_str());
		return;
	}
}
	
//---------------------------------------------------------------------------------------------------------------

//opens an ASC file for writing and places the header records

FILE * openasc_w(std::string filename, asc_header hdr)
{
	FILE * ofile=fopen(filename.c_str(),"w");
	if(!ofile){
		std::cout<<"Error: Failed to open output file (\""<<filename.c_str()<<"\")."<<std::endl;
		return NULL;
	}
	
	fprintf(ofile,"ncols %d\n",hdr.ncols);
	fprintf(ofile,"nrows %d\n",hdr.nrows);
	fprintf(ofile,"xllcorner %lf\n",hdr.xll);
	fprintf(ofile,"yllcorner %lf\n",hdr.yll);
	fprintf(ofile,"cellsize %lf\n",hdr.cellsize);
	fprintf(ofile,"nodata_value %lf\n",hdr.nodata_value);
	
	return ofile;
}

//---------------------------------------------------------------------------------------------------------------

//opens an ASC file for writing and places the header records

FILE * openasc_r(std::string filename, asc_header &hdr)
{
	FILE * ifile=fopen(filename.c_str(),"r");
	if(!ifile){
		std::cout<<"Error: Failed to open \""<<filename.c_str()<<"\"."<<std::endl;
		return NULL;
	}
	
	//READ HEADER
	char * cbuf=new char [50];
	fscanf(ifile,"%s",cbuf);	
		
	//check file type
	if(strncasecmp(cbuf,"ncols",5)!=0){
		fclose(ifile);
		delete [] cbuf;
		std::cout<<"Error: File \""<<filename.c_str()<<"\" is not a valid Esri (tm) ASCII grid file."<<std::endl;
		return NULL;
	}
		
	fscanf(ifile,"%i",&(hdr.ncols));
	fscanf(ifile,"%s",cbuf);
	fscanf(ifile,"%i",&(hdr.nrows));
	fscanf(ifile,"%s",cbuf);
	fscanf(ifile,"%lf",&(hdr.xll));
	fscanf(ifile,"%s",cbuf);
	fscanf(ifile,"%lf",&(hdr.yll));
	fscanf(ifile,"%s",cbuf);
	fscanf(ifile,"%lf",&(hdr.cellsize));
	fscanf(ifile,"%s",cbuf);
	fscanf(ifile,"%lf",&(hdr.nodata_value));
	
	return ifile;
}

//---------------------------------------------------------------------------------------------------------------

point makepoint(double x, double y)
{ 
	point p; 
	p.x=x; 
	p.y=y; 
	return p; 
}
 