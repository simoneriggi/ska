// ***********************************************************************
// * License and Disclaimer                                              *
// *                                                                     *
// * Copyright 2016 Simone Riggi																			   *
// *																																	   *
// * This file is part of Caesar. 																		   *
// * Caesar is free software: you can redistribute it and/or modify it   *
// * under the terms of the GNU General Public License as published by   *
// * the Free Software Foundation, either * version 3 of the License,    *
// * or (at your option) any later version.                              *
// * Caesar is distributed in the hope that it will be useful, but 			 *
// * WITHOUT ANY WARRANTY; without even the implied warranty of          * 
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                *
// * See the GNU General Public License for more details. You should     * 
// * have received a copy of the GNU General Public License along with   * 
// * Caesar. If not, see http://www.gnu.org/licenses/.                   *
// ***********************************************************************
/**
* @file GraphicsUtils.h
* @class GraphicsUtils
* @brief Utility functions for graphics tasks
*
* Utility functions for graphics tasks
* @author S. Riggi
* @date 15/01/2016
*/


#ifndef GraphicsUtils_h
#define GraphicsUtils_h 1

#include <TObject.h>
#include <TGaxis.h>
#include <TPolyLine.h>

namespace Caesar {

class Img;


enum ColorPaletteStyle {
	eRAINBOW= 0,
	eBLACKWHITE= 1,
	eBLACKBODY= 2,
	eHOT2COLD= 3,
	eCOLD2HOT= 4,
	eTHERMAL= 5
};

class GraphicsUtils : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    GraphicsUtils();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~GraphicsUtils();

		
	public:
		static int SetThermalPalette(int ncolors=999);
		static int SetHotColdPalette(int ncolors=999);
		static int SetColdHotPalette(int ncolors=999);
		static int SetBWPalette(int ncolors=999);
		static int UpdateGAxis();
		static int SetWCSAxis(Img* img,TGaxis& xaxis,TGaxis& yaxis,int coordSystem=-1);
		static int SetWCSProjGrid(Img* img,std::vector<TPolyLine>& gridx,std::vector<TPolyLine>& gridy,int coordSystem);
		static Img* FindImageFromPad();
		static int PadUpdater();

	private:
	
		ClassDef(GraphicsUtils,1)
};

#ifdef __MAKECINT__
#pragma link C++ class GraphicsUtils+;
//#pragma link C++ enum ColorPaletteStyle+;
#endif	

}//close namespace


#endif 
