/*******************************************************************************
 * Copyright (c) 2014 Chris Satterwhite
 * Copyright (c) 2018 David D. Marshall <ddmarsha@calpoly.edu>
 *
 * This program and the accompanying materials are made
 * available under the terms of the Eclipse Public License 2.0
 * which is available at https://www.eclipse.org/legal/epl-2.0/
 *
 * See LICENSE.md file in the project root for full license information.
 *
 * SPDX-License-Identifier: EPL-2.0
 *
 * Contributors:
 *    Chris Satterwhite - initial code and implementation
 *    Connor Sousa - Vortex particle implementation
 *    David D. Marshall - misc. changes
 ******************************************************************************/

#include "geometry.h"
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#define _USE_MATH_DEFINES
#include <math.h>

// Destructor //

geometry::~geometry()
{

    for (surfaces_index_type i=0; i<surfaces.size(); i++)
    {
        delete surfaces[i];
    }
    surfaces.clear();
    for (wakes_index_type i=0; i<wakes.size(); i++)
    {
        delete wakes[i];
    }
    wakes.clear();
    for (edges_index_type i=0; i<edges.size(); i++)
    {
        delete edges[i];
    }
    edges.clear();
    for (nodes_index_type i=0; i<nodes.size(); i++)
    {
        delete nodes[i];
    }
    nodes.clear();
    
}

// Copy Constructor //

geometry::geometry(const geometry& copy)
  : pOctree(copy.pOctree), nNodes(copy.nNodes), nTris(copy.nTris), A(copy.A), B(copy.B), C(copy.C),
	writeCoeffFlag(false), vortPartFlag(false), infCoeffFile(copy.infCoeffFile), dt(copy.dt),
	inputV(copy.inputV)
{
    for (surfaces_index_type i=0; i<copy.surfaces.size(); i++)
    {
        surfaces[i] = new surface(*copy.surfaces[i]);
    }
    for (wakes_index_type i=0; i<copy.wakes.size(); i++)
    {
        wakes[i] = new wake(*copy.wakes[i]);
    }
    for (bodyPanels_index_type i=0; i<copy.bPanels.size(); i++)
    {
        bPanels[i] = new bodyPanel(*copy.bPanels[i]);
    }
    for (wakePanels_index_type i=0; i<copy.wPanels.size(); i++)
    {
        wPanels[i] = new wakePanel(*copy.wPanels[i]);
    }
    for (nodes_index_type i=0; i<copy.nodes.size(); i++)
    {
        nodes[i] = new cpNode(*copy.nodes[i]);
    }
    for (edges_index_type i=0; i<copy.edges.size(); i++)
    {
        edges[i] = new edge(*copy.edges[i]);
    }
}

// Assignment Operator //

geometry& geometry::operator=(const geometry &rhs)
{
    if (this == &rhs)
    {
        return (*this);
    }

    pOctree = rhs.pOctree;
    nNodes = rhs.nNodes;
    nTris = rhs.nTris;
    A = rhs.A;
    B = rhs.B;
    infCoeffFile = rhs.infCoeffFile;
    

    for (surfaces_index_type i=0; i<rhs.surfaces.size(); i++)
    {
        surfaces[i] = new surface(*rhs.surfaces[i]);
    }
    for (wakes_index_type i=0; i<rhs.wakes.size(); i++)
    {
        wakes[i] = new wake(*rhs.wakes[i]);
    }
    for (bodyPanels_index_type i=0; i<rhs.bPanels.size(); i++)
    {
        bPanels[i] = new bodyPanel(*rhs.bPanels[i]);
    }
    for (wakePanels_index_type i=0; i<rhs.wPanels.size(); i++)
    {
        wPanels[i] = new wakePanel(*rhs.wPanels[i]);
    }
    for (nodes_index_type i=0; i<rhs.nodes.size(); i++)
    {
        nodes[i] = new cpNode(*rhs.nodes[i]);
    }
    for (edges_index_type i=0; i<rhs.edges.size(); i++)
    {
        edges[i] = new edge(*rhs.edges[i]);
    }

    return *this;
}



void geometry::supSortSurfs(std::string tri_file)
{
	// open tri file
	std::ifstream fid;
	fid.open(tri_file);
	if (fid.is_open())
	{
		// read and assign number of nodes and number of tris
		fid >> nNodes >> nTris;
		nNodesOrig = nNodes;

		// compute panel normals

		// Read XYZ Locations of Nodes

		//std::vector<Eigen::Vector3d> tempPnts;
		//Eigen::Vector3d tempPnt;

		//struct tempPntsStruct
		//{
		//	Eigen::Vector3d tempPnt;
		//	size_t index;
		//	bool isCopied = false;
		//};

		std::vector<tempPntsStruct> tempPnts;

		for (size_t i = 0; i<nNodes; i++)
		{
			tempPnts.push_back(tempPntsStruct());
			fid >> tempPnts[i].tempPnt(0) >> tempPnts[i].tempPnt(1) >> tempPnts[i].tempPnt(2);
			//tempPnts[i].index = i;

			//fid >> tempPnt(0) >> tempPnt(1) >> tempPnt(2);
			//tempPnts.push_back(tempPnt);
		}

		// read connectivity
		
		std::vector<Eigen::Vector3d> tempCons;
		Eigen::Vector3d tempCon;

		/*struct tempCons
		{
			std::vector<Eigen::Vector3d> tempCons;
			bool isCopied = false;
		};*/



		for (size_t i = 0; i<nTris; i++)
		{
			fid >> tempCon(0) >> tempCon(1) >> tempCon(2);
			//fid >> tempCon(0) >> tempCon(2) >> tempCon(1);
			tempCons.push_back(tempCon);
			tempCons[i].x() -= 1;
			tempCons[i].y() -= 1;
			tempCons[i].z() -= 1;
		}

		// compute normals
		std::vector<Eigen::Vector3d> tempNorms;
		//Eigen::Vector3d tempNorm;
		Eigen::Vector3d p0, p1, p2, a, b;
		std::vector<size_t> upperPans, lowerPans, otherPans;
		for (size_t i = 0; i < nTris; i++)
		{
			p0 = tempPnts[tempCons[i].x()].tempPnt;
			p1 = tempPnts[tempCons[i].y()].tempPnt;
			p2 = tempPnts[tempCons[i].z()].tempPnt;

			/*p0 = tempPnts[tempCons[i].x()];
			p1 = tempPnts[tempCons[i].y()];
			p2 = tempPnts[tempCons[i].z()];*/

			a = p1 - p0;
			b = p2 - p0;

			tempNorms.push_back(a.cross(b));
			tempNorms[i].normalize();

			///////////////////////////////////////////////////////////////////////////
			//std::cout << tempNorms[i] << "\n" << std::endl;
			///////////////////////////////////////////////////////////////////////////

			// check sign of normal z-component (i.e. upper or lower panel)
			if (tempNorms[i].z() > 1.0e-10) // upper
			{
				upperPans.push_back(i);
			}
			else if (tempNorms[i].z() < -1.0e-10) // lower
			{
				lowerPans.push_back(i);
			}
			else // other
			{
				otherPans.push_back(i);
			}
		}

		size_t nUpper, nLower, nOther;
		nUpper = upperPans.size();
		nLower = lowerPans.size();
		nOther = otherPans.size();

		//std::vector<size_t> origCopiedPnts, copiedPnts;

		//std::vector<int> commonCons;

		// Compare UPPERPANS w/ LOWERPANS and find common points, make copies and assign
		supSort(upperPans, lowerPans, tempPnts, tempCons, nNodes);

		// Saturday night: Multiple copies made of same point. Need to figure out how to make just one copy of needed points
		//				   and apply them as needed

		// Create a struct for each point which has node and a flag which says whether or not it's been copied already

		////////	 // NOT WORKING
		//		ANSWER?? : Make two vectors. One for original points, and one for copies of original points. Make
		//				   indices align between two vectors. Then check vectors each time through



		// ARE THESE NEEDED  ?????????????????????????????????????????

		/*
		// Compare UPPERPANS w/ OTHERPANS " "
		supSort(upperPans, otherPans, tempPnts, tempCons, nNodes);

		// Compare LOWERPANS w/ OTHERPANS " "
		supSort(lowerPans, otherPans, tempPnts, tempCons, nNodes);
		*/



		// read panel IDs (so that they can be written back into the new .tri file)
		int tempID;
		std::vector<int> tempIDs;
		for (size_t i = 0; i < nTris; i++)
		{
			fid >> tempID;
			tempIDs.push_back(tempID);
		}
		fid.close();


		// open and clear new .tri file
		fid.open(tri_file, std::ofstream::out | std::ofstream::trunc);
		fid.close();

		// write new panel and node data to new .tri file
		std::ofstream fid;
		fid.open(tri_file);

		// write nNodes and nTris
		fid << nNodes << " " << nTris << "\n";

		// write nodes
		for (size_t i = 0; i < nNodes; i++)
		{
			fid << "\t" << tempPnts[i].tempPnt.x() << "\t" << tempPnts[i].tempPnt.y() << "\t" << tempPnts[i].tempPnt.z() << "\n";
		}

		// write connectivity
		for (size_t i = 0; i < nTris; i++)
		{
			fid << tempCons[i].x()+1 << " " << tempCons[i].y()+1 << " " << tempCons[i].z()+1 << "\n";
		}

		// write panel IDs
		for (size_t i = 0; i < nTris; i++)
		{
			fid << tempIDs[i] << "\n";
		}


		fid.close();

	}
	else
	{
		std::cout << "ERROR : Geometry file not found" << std::endl;
		exit(EXIT_FAILURE);
	}
}


//void geometry::supSort(std::vector<size_t> &panSet1, std::vector<size_t> &panSet2, std::vector<Eigen::Vector3d> &tempPnts, std::vector<Eigen::Vector3d> &tempCons, size_t &nNodes, std::vector<int> &commonCons)
void geometry::supSort(std::vector<size_t> &panSet1, std::vector<size_t> &panSet2, std::vector<tempPntsStruct> &tempPnts, std::vector<Eigen::Vector3d> &tempCons, size_t &nNodes)
{

	size_t nPanSet1 = panSet1.size();
	size_t nPanSet2 = panSet2.size();

	size_t count = 0;

	for (size_t i = 0; i < nPanSet1; i++)
	{
		for (size_t j = 0; j < nPanSet2; j++)
		{
			//std::vector<int> commonCon;
			std::vector<int> commonCons;

			// put panel sets 1 and 2 connectivity data into STL vectors
			std::vector<double> compCons1 = { tempCons[panSet1[i]].x(), tempCons[panSet1[i]].y(), tempCons[panSet1[i]].z() };
			std::vector<double> compCons2 = { tempCons[panSet2[j]].x(), tempCons[panSet2[j]].y(), tempCons[panSet2[j]].z() };

			// sort the two vectors
			std::sort(compCons1.begin(), compCons1.end());
			std::sort(compCons2.begin(), compCons2.end());

			// initialise a vector store the common values and an iterator to traverse this vector 
			std::vector<int> v(compCons1.size() + compCons2.size());
			std::vector<int>::iterator it, st;

			// define iterator
			it = std::set_intersection(compCons1.begin(), compCons1.end(), compCons2.begin(), compCons2.end(), v.begin());

			// find common elements
			for (st = v.begin(); st != it; ++st)
			{
				//commonCon.push_back(*st);
				//commonCon = *st;
				commonCons.push_back(*st);
			}

			// check if there were any common elements
			if (commonCons.size() > 0)
			//if (commonCon.size() > 0)
			{
				for (size_t k = 0; k < commonCons.size(); k++)
				{
					// check if common point has already been copied

					// from common connectivity index, go into points and check isCopied flag
					if (tempPnts[commonCons[k]].isCopied == true)
					{
						if (tempCons[panSet2[j]].x() == commonCons[k])
						{
							tempCons[panSet2[j]].x() = tempPnts[commonCons[k]].copyIndex;
							//tempCons[panSet2[j]].x() = nNodes - 1;
						}
						else if (tempCons[panSet2[j]].y() == commonCons[k])
						{
							tempCons[panSet2[j]].y() = tempPnts[commonCons[k]].copyIndex;
						}
						else if (tempCons[panSet2[j]].z() == commonCons[k])
						{
							tempCons[panSet2[j]].z() = tempPnts[commonCons[k]].copyIndex;
						}
					}
					else
					{
						nNodes += 1;

						// copy common points to the end of the list of points
						tempPnts.push_back(tempPntsStruct());
						tempPnts[nNodes-1].tempPnt = tempPnts[commonCons[k]].tempPnt;
						tempPnts[commonCons[k]].isCopied = true;
						tempPnts[commonCons[k]].copyIndex = nNodes - 1;

						// replace index of common point in the LOWER panel with index of new copy of common point
						if (tempCons[panSet2[j]].x() == commonCons[k])
						{
							tempCons[panSet2[j]].x() = nNodes - 1;	// same as copyIndex, doesn't matter which is used
						}
						if (tempCons[panSet2[j]].y() == commonCons[k])
						{
							tempCons[panSet2[j]].y() = nNodes - 1;
						}
						if (tempCons[panSet2[j]].z() == commonCons[k])
						{
							tempCons[panSet2[j]].z() = nNodes - 1;
						}
					}
					


					/*size_t m = 0;
					while(commonCon[k] !=)

					for (size_t m = 0; m < commonCons.size(); m++)
					{

					}*/

					//// copy common points to the end of the list of points
					//tempPnts.push_back(tempPnts[commonCons[k]]);

					//// replace index of common point in the LOWER panel with index of new copy of common point
					//if (tempCons[panSet2[j]].x() == commonCons[k])
					//{
					//	nNodes += 1;
					//	tempCons[panSet2[j]].x() = nNodes - 1;
					//}
					//if (tempCons[panSet2[j]].y() == commonCons[k])
					//{
					//	nNodes += 1;
					//	tempCons[panSet2[j]].y() = nNodes - 1;
					//}
					//if (tempCons[panSet2[j]].z() == commonCons[k])
					//{
					//	nNodes += 1;
					//	tempCons[panSet2[j]].z() = nNodes - 1;
					//}
				}
			}

		}
	}


}



// upperPans VS lowerPans

//for (size_t i = 0; i < nUpper; i++)
//{
//	for (size_t j = 0; j < nLower; j++)
//	{
//		std::vector<int> commonCons;

//		// put upper and lower panel connectivity data into STL vectors
//		std::vector<double> compCons1 = { tempCons[upperPans[i]].x(), tempCons[upperPans[i]].y(), tempCons[upperPans[i]].z() };
//		std::vector<double> compCons2 = { tempCons[lowerPans[j]].x(), tempCons[lowerPans[j]].y(), tempCons[lowerPans[j]].z() };

//		/*std::vector<double> compCons1_Orig = { tempCons[upperPans[i]].x(), tempCons[upperPans[i]].y(), tempCons[upperPans[i]].z() };
//		std::vector<double> compCons2_Orig = { tempCons[lowerPans[j]].x(), tempCons[lowerPans[j]].y(), tempCons[lowerPans[j]].z() };
//		std::vector<double> compCons1 = compCons1_Orig;
//		std::vector<double> compCons2 = compCons2_Orig;*/

//		// sort the two vectors
//		std::sort(compCons1.begin(), compCons1.end());
//		std::sort(compCons2.begin(), compCons2.end());

//		// initialise a vector store the common values and an iterator to traverse this vector 
//		std::vector<int> v(compCons1.size() + compCons2.size());
//		std::vector<int>::iterator it, st;

//		// define iterator
//		it = std::set_intersection(compCons1.begin(), compCons1.end(), compCons2.begin(), compCons2.end(), v.begin());

//		// find common elements
//		for (st = v.begin(); st != it; ++st)
//		{
//			commonCons.push_back(*st);
//		}

//		// check if there were any common elements
//		if (commonCons.size() > 0)
//		{
//			for (size_t k = 0; k < commonCons.size(); k++)
//			{
//				// copy common points to the end of the list of points
//				tempPnts.push_back(tempPnts[commonCons[k]]);

//				// replace index of common point in the LOWER panel with index of new copy of common point
//				if (tempCons[lowerPans[j]].x() == commonCons[k])
//				{
//					nNodes += 1;
//					tempCons[lowerPans[j]].x() = nNodes - 1;
//				}
//				if (tempCons[lowerPans[j]].y() == commonCons[k])
//				{
//					nNodes += 1;
//					tempCons[lowerPans[j]].y() = nNodes - 1;
//				}
//				if (tempCons[lowerPans[j]].z() == commonCons[k])
//				{
//					nNodes += 1;
//					tempCons[lowerPans[j]].z() = nNodes - 1;
//				}
//			}
//		}

//	}
//}



void geometry::readTri(std::string tri_file, bool normFlag)
{
    std::ifstream fid;
    fid.open(tri_file);
    if (fid.is_open())
    {
		// Check if Mach number is in the transonic regime and ask for user input
		if (0.6 < inputMach && inputMach < 1.3)
		{
			std::string in;
			std::cout << "WARNING: Input Mach number could result in transonic flow phenomena which will not be modeled properly!" << std::endl;
			std::cout << "\tWould you like to proceed anyway?" << std::endl;
			std::cout << "\t\t< Y > - Yes, proceed any." << std::endl;
			std::cout << "\t\t< N > - No, exit program." << std::endl;
			std::cin >> in;
			std::cout << std::endl;
			if (in == "N" || in == "n")
			{
				exit(EXIT_FAILURE);
			}
		}



		//---------------------------------------------------------------------------------------------------------------------//
		//---------------------------------------------------------------------------------------------------------------------//
		//---------------------------------------------------------------------------------------------------------------------//
		//---------------------------------------------------------------------------------------------------------------------//

		// NOTES: - number of panels remains the same, it's just the node count that changes, and where the nodes are assigned

		// Algorithm for solving supersonic solver leading and trailing edge problems
		//		1. Read .tri file
		//		2. Identify leading and trailing edge panels/nodes
		//		3. Make copies of leading/trailing edge nodes
		//		4. Assign one copy of node to upper panel, and the other copy to the lower panel
		//		5. I think that's it??

		//// Create temporary fid
		//std::ifstream fidTemp;
		//fidTemp.open(tri_file);

		//// Read first line of tri file
		//size_t nNodesTemp, nTrisTemp;
		//fidTemp >> nNodes >> nTris;

		//// Read XYZ Locations of Nodes and Create temporary nodes
		//Eigen::Vector3d pntTemp;
		//nodes_type nodesTemp;
		//cpNode* n;
		//for (size_t i = 0; i<nNodes; i++)
		//{
		//	fid >> pntTemp(0) >> pntTemp(1) >> pntTemp(2);
		//	n = new cpNode(pntTemp, i);
		//	nodesTemp.push_back(n);
		//}


		// Need to only go through body nodes, and not wake nodes --> use IDs

		//---------------------------------------------------------------------------------------------------------------------//
		//---------------------------------------------------------------------------------------------------------------------//
		//---------------------------------------------------------------------------------------------------------------------//
		//---------------------------------------------------------------------------------------------------------------------//



        std::cout << "Reading Geometry..." << std::endl;
        fid >> nNodes >> nTris;
        Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> connectivity(nTris,3);
        Eigen::VectorXi allID(nTris);
        std::vector<int> surfIDs;
        std::vector<int> wakeIDs;
        std::vector<int> surfTypes;


        // Read XYZ Locations of Nodes
        Eigen::Vector3d pnt;
        cpNode* n;
        for (size_t i=0; i<nNodes; i++)
        {
            fid >> pnt(0) >> pnt(1) >> pnt(2);
            n = new cpNode(pnt,i);
            nodes.push_back(n);
        }


        // Temporarily Store Connectivity
		if (inputMach > 1.0)
		{
			// For supersonic scheme, need to go around panel CCW -> switch 1 and 2 below relative to original
			for (Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index i = 0; i<static_cast<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index>(nTris); i++)
			{
				fid >> connectivity(i, 0) >> connectivity(i, 2) >> connectivity(i, 1);
			}
		}
		else
		{
			for (Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index i = 0; i<static_cast<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index>(nTris); i++)
			{
				fid >> connectivity(i, 0) >> connectivity(i, 1) >> connectivity(i, 2);
			}
		}

        connectivity = connectivity.array()-1; //Adjust for 0 based indexing

        // Scan Surface IDs and collect Unique IDs
        size_t wakeNodeStart = nNodes;
        size_t wakeTriStart = nTris;
        for (Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index i=0; i<static_cast<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index>(nTris); i++)
        {
            fid >> allID(i);
            if (i == 0 || allID(i) != allID(i-1))
            {
                if (allID(i) > 1000)
                {
                    wakeIDs.push_back(allID(i));
                }
                else
                {
                    surfIDs.push_back(allID(i));
                }
            }
            if (allID(i) > 1000 && allID(i-1) < 1000)
            {
                wakeNodeStart = connectivity.row(i).minCoeff();
                wakeTriStart = static_cast<size_t>(i);
            }
        }

        if (wakeIDs.size() > 0)
        {
            correctWakeConnectivity(wakeNodeStart, wakeTriStart, connectivity);
        }

        // Read in Normals if included in input file
        Eigen::MatrixXd norms = Eigen::MatrixXd::Zero(static_cast<Eigen::MatrixXd::Index>(nTris),3);
        if (normFlag)
        {
            std::cout << "Reading Bezier Normals from Geometry File..." << std::endl;
            for (Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index i=0; i<static_cast<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index>(nTris); i++)
            {
                fid >> norms(i,0) >> norms(i,1) >> norms(i,2);
            }
        }

        std::cout << "Generating Panel Geometry..." << std::endl;
        
        createSurfaces(connectivity,norms,allID);
        
        /* Addition of 2 row buffer wake
         - Calculate time step (for panel length) (only if not included)
         - Find and delete any included wake surfaces
         - Create buffer wake
         */
        
        if(vortPartFlag)
        {
            // VSP tags wakes with surface ID starting at 1000
            for(int i=0; i<allID.size(); i++)
            {
                if(allID(i) >= 1000)
                {
                    std::cout << "ERROR: Please use geometry file without wake panels when using the vortex particle wake option." << std::endl;
                    std::exit(0);
                }
            }
            
            // Calculate timestep
            calcTimeStep();
            
            // Find existing wake surfaces and delete
            while ( wakes.size() > 0 ) {
                wakes.erase(wakes.begin(), wakes.end());
            }
            
            // Build wake panel nodes:
            
            // Finding trailing edges and nodes
            std::vector< edge* > TEedges;
            std::vector< cpNode* > TEnodes;
            
            for(edges_index_type i=0; i<edges.size(); i++ ){
                if( edges[i]->isTE() ){
                    TEedges.push_back(edges[i]);
                }
            }
            for(nodes_index_type i=0; i<nodes.size(); i++ ){
                if( nodes[i]->isTE() ){
                    TEnodes.push_back(nodes[i]);
                }
            }
            
            Eigen::MatrixXd wakeNodes(2*TEnodes.size(),3);
            Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> wakeConnectivity(2*TEedges.size(),4);
            std::vector<size_t> VPwakeID, newNodesIndex, usedTENodesIndex;
            size_t nodeCounter=0, panelCounter=0;
            
            for(size_t i=0; i<TEedges.size(); i++)
            {
                
                // Find the trialing nodes (by index?)
                size_t n1index = TEedges[i]->getN1()->getIndex();
                size_t n2index = TEedges[i]->getN2()->getIndex();
                size_t n1firstIndex, n1secIndex, n2firstIndex, n2secIndex;
                
                bool isUsed = false;
                for(size_t j=0; j<usedTENodesIndex.size(); j++)
                {
                    if(n1index == usedTENodesIndex[j])
                    {
                        isUsed = true;
                        n1firstIndex = newNodesIndex[2*j]+nNodes;
                        n1secIndex = newNodesIndex[2*j+1]+nNodes;
                    }
                }
                
                if(isUsed == false)
                { // If it's used, don't change the n1 index and then get the first and sec. node indices.
                    usedTENodesIndex.push_back(n1index);
                    
                    double VinfLocal = inputV;
                    
                    wakeNodes.row(static_cast<Eigen::MatrixXd::Index>(nodeCounter)) = nodes[n1index]->firstProjNode(dt, VinfLocal);
                    newNodesIndex.push_back(nodeCounter);
                    n1firstIndex = nodeCounter+nNodes;
                    nodeCounter++;
                    
                    wakeNodes.row(static_cast<Eigen::MatrixXd::Index>(nodeCounter)) = nodes[n1index]->secProjNode(dt, VinfLocal);
                    newNodesIndex.push_back(nodeCounter);
                    n1secIndex = nodeCounter+nNodes;
                    nodeCounter++;
                }
                //N2
                isUsed = false;
                for(size_t j=0; j<usedTENodesIndex.size();j++)
                {
                    if(n2index == usedTENodesIndex[j]){
                        isUsed = true;
                        n2firstIndex = newNodesIndex[2*j]+nNodes;
                        n2secIndex = newNodesIndex[2*j+1]+nNodes;
                    }
                }
                
                if(isUsed == false)
                {
                    usedTENodesIndex.push_back(n2index);
                    
                    double VinfLocal = inputV;
                    wakeNodes.row(static_cast<Eigen::MatrixXd::Index>(nodeCounter)) = nodes[n2index]->firstProjNode(dt, VinfLocal);
                    newNodesIndex.push_back(nodeCounter);
                    n2firstIndex = nodeCounter+nNodes;
                    nodeCounter++;
                    
                    wakeNodes.row(static_cast<Eigen::MatrixXd::Index>(nodeCounter)) = nodes[n2index]->secProjNode(dt, VinfLocal);
                    newNodesIndex.push_back(nodeCounter);
                    n2secIndex = nodeCounter+nNodes;
                    nodeCounter++;
                }
                
                
                wakeConnectivity.row(static_cast<Eigen::MatrixXd::Index>(panelCounter)) << n1index, n2index, n2firstIndex, n1firstIndex; //built TE first
                
                panelCounter++;
                VPwakeID.push_back(1001); // First row of wake panels
                isFirstPanel.push_back(true);
                
                wakeConnectivity.row(static_cast<Eigen::MatrixXd::Index>(panelCounter)) << n1firstIndex, n2firstIndex, n2secIndex, n1secIndex;
                panelCounter++;
                VPwakeID.push_back(1001); // Second row
                isFirstPanel.push_back(false);
            }
            
            
            // Append nodes and adjust connectivity
            for (Eigen::MatrixXd::Index i=0; i<wakeNodes.rows(); i++)
            {
                n = new cpNode(wakeNodes.row(i),static_cast<size_t>(i)+nNodes);
                nodes.push_back(n);
            }
            Eigen::MatrixXd wakeNorms = Eigen::MatrixXd::Zero(wakeConnectivity.rows(),3);
            createVPWakeSurfaces(wakeConnectivity,wakeNorms,VPwakeID,isFirstPanel);
            
            nNodes = nodes.size();
            nTris += wakes[0]->getPanels().size(); // Include buffer wake
            
        }
        
        
        std::cout << "\tNodes : " << nodes.size() << std::endl;
        std::cout << "\tEdges : " << edges.size() << std::endl;
        std::cout << "\tPanels : " << nTris << std::endl;

        // Erase duplicate node pointers
        std::sort( nodes.begin(), nodes.end() );
        nodes.erase( std::unique( nodes.begin(), nodes.end() ), nodes.end() );

        for (nodes_index_type i=0; i<nodes.size(); i++)
        {
            nodes[i]->setIndex(i);
        }

        std::cout << "Building Octree..." << std::endl;
        
        createOctree();

        // Set neighbors
        std::cout << "Finding Panel Neighbors..." << std::endl;

        for (edges_index_type i=0; i<edges.size(); i++)
        {
            edges[i]->setNeighbors();
        }


        bool wakeMergeFlag = false;
        if (wakes.size() > 1)
        {
            std::vector<wake*> newWakes;
            for (wakes_index_type i=0; i<wakes.size(); i++)
            {
                for (wakes_index_type j=i; j<wakes.size(); j++)
                {
                    if (wakes[i]->isSameWake(wakes[j]))
                    {
                        wakes[i]->mergeWake(wakes[j]);
                        delete wakes[j];
                        newWakes.push_back(wakes[i]);
                        wakeMergeFlag = true;
                    }
                }
            }
            if (wakeMergeFlag == true)
            {
                wakes = newWakes;
            }
        }


        // Collect all panels in geometry

        std::vector<bodyPanel*> tempB;
        std::vector<wakePanel*> tempW;
        for (surfaces_index_type i=0; i<surfaces.size(); i++)
        {
            tempB = surfaces[i]->getPanels();
            bPanels.insert(bPanels.begin(),tempB.begin(),tempB.end());
        }
        for (wakes_index_type i=0; i<wakes.size(); i++)
        {
            tempW = wakes[i]->getPanels();
            for (size_t j=0; j<tempW.size(); j++)
            {
                if (tempW[j]->isSecondRow == false)
                {
                    wPanels.push_back(tempW[j]);
                }
            }
        }

        
		// Only need tip patch identification for lower-order scheme
		if (!subHOMFlag && inputMach < 1.0)
		{
			// Check panels for tip patches.  Needed to do 2D CHTLS to avoid nonphysical results near discontinuity at trailing edge.
			for (bodyPanels_index_type i = 0; i<bPanels.size(); i++)
			{
				bPanels[i]->setTipFlag();
			}
			for (bodyPanels_index_type i = 0; i<bPanels.size(); i++)
			{
				bPanels[i]->setCluster();
			}
		}


		if (inputMach > 1.0)
		{
			// For supersonic scheme, flip panel normal due to change in ordering of panel vertices

			for (bodyPanels_index_type i = 0; i<bPanels.size(); i++)
			{
				bPanels[i]->supFlipNormal();
			}

			// Check for superinclined panels. Ask if user wants to continue if at least one is found

			std::cout << "Checking for superinclined panels...";
			bool isSupInclined = false;
			Eigen::Vector3d nWind;
			double B = sqrt(pow(inputMach, 2) - 1);
			setBodyToWind(alpha, beta);

			size_t count = 0;
			for (size_t i = 0; i < bPanels.size(); i++)
			{
				isSupInclined = bPanels[i]->supSuperinclinedCheck(B, body2wind); // returns true if sup-inclined panel is found
				if (isSupInclined)
				{
					count += 1;
				}
			}
			if (count > 0)
			{
				std::string in;
				std::cout << "\n\n" << count << " superinclined panels found, input geometry must be modified or results will be in error." << std::endl;
				std::cout << "\tWould you like to proceed anyway?" << std::endl;
				std::cout << "\t\t< Y > - Yes, proceed any." << std::endl;
				std::cout << "\t\t< N > - No, exit program." << std::endl;
				std::cin >> in;
				std::cout << std::endl;
				if (in == "N" || in == "n")
				{
					exit(EXIT_FAILURE);
				}
			}
			else
			{
				std::cout << " None found" << std::endl;
			}
		}


		if (subHOMFlag || inputMach > 1.0)
		{
			// Get control point data for linear schemes
			Eigen::Vector3d windDir = supComputeWindDir();
			for (nodes_index_type i = 0; i < nodes.size(); i++)
			{
				if (nodes[i]->getBodyPans().size() > 0)
				{
					bodyNodes.push_back(nodes[i]);
					nodes[i]->setLinCPnormal();

					if (i > nNodesOrig - 1)
					{
						nodes[i]->setLinCPoffsetCopy();
					}
					else
					{
						nodes[i]->setLinCPoffset();
					}
					//nodes[i]->setLinCPoffset();
				}
				else
				{
					wakeNodes.push_back(nodes[i]);
				}
			}


			//----------------------------------------------------------------------------------------------------//

			// Organize body nodes for wake handling... wake handling is postponed for the time being 10/18/2018

			/*nodes_index_type j = 0;
			nodes_index_type k = bodyNodes.size();
			nodes_type tempNodes;
			tempNodes.resize(static_cast<Eigen::MatrixXd::Index>(nodes.size()));
			for (nodes_index_type i = 0; i < nodes.size(); i++)
			{
				if (nodes[i]->getBodyPans().size() > 0)
				{
					tempNodes[j] = nodes[i];
					tempNodes[j]->setLinCPnormal();
					tempNodes[j]->setLinCPoffset();
					j += 1;
				}
				else
				{
					tempNodes[k] = nodes[i];
					k += 1;
				}
			}
			nodes = tempNodes;

			for (size_t i = 0; i < nodes.size(); i++)
			{
				nodes[i]->setIndex(static_cast<int>(i));
			}*/

			//----------------------------------------------------------------------------------------------------//
		}


        // Calculate influence coefficient matrices

        bool read = false;

        if (infCoeffFileExists())
        {
            std::string in;
            std::cout << "\nInfluence Coefficients have already been calculated for a geometry with this name, would you like to use these coefficients?" << std::endl;
            std::cout << "\t< Y > - Yes, use coefficients." << std::endl;
            std::cout << "\t< N > - No, recalculate them." << std::endl;
            std::cin >> in;
            std::cout << std::endl;
            if (in == "Y" || in == "y")
            {
                readInfCoeff();
                read = true;
                std::cout << "take this out" << std::endl;
            }
        }

        if (!read)
        {
            std::cout << "Building Influence Coefficient Matrix..." << std::endl;
			if (subHOMFlag) // linear doublet subsonic scheme
			{
				linSetInfCoeff();
			}
			else if (inputMach > 1.0) // supersonic scheme
			{
				supSetInfCoeff();
			}
			else // original subsonic constant doublet scheme
			{
				setInfCoeff();
			}
        }
    }
    else
    {
        std::cout << "ERROR : Geometry file not found" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void geometry::correctWakeConnectivity(size_t wakeNodeStart,size_t wakeTriStart,Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> &connectivity)
{
    Eigen::Vector3d vec;
    Eigen::Matrix<size_t,Eigen::Dynamic,2> indReps; // [toReplace, replaceWith]
    double tol = shortestEdge(connectivity);
    int count = 0;
    for (size_t i=0; i<wakeNodeStart; i++)
    {
        for (size_t j=wakeNodeStart; j<nNodes; j++)
        {
            vec = nodes[i]->getPnt()-nodes[j]->getPnt();
            if (vec.lpNorm<Eigen::Infinity>() < tol)
            {

                nodes[j] = nodes[i];
                count++;
                indReps.conservativeResize(count,2);
                indReps(count-1,0) = j;
                indReps(count-1,1) = i;
            }
        }
    }

    for (Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index i=static_cast<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index>(wakeTriStart); i<static_cast<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index>(nTris); i++)
    {
        for (int j=0; j<connectivity.cols(); j++)
        {
            for (int k=0; k<indReps.rows(); k++)
            {
                if (connectivity(i,j) == indReps(k,0))
                {
                    connectivity(i,j) = indReps(k,1);
                }
            }
        }
    }
}

double geometry::shortestEdge(const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> &connectivity)
{
    int nrows,ncols;
    Eigen::Vector3d p1,p2;
    nrows = (int)connectivity.rows();
    ncols = (int)connectivity.cols();
    double l;
    double shortest = 100000;
    for (int i=0; i<nrows; i++)
    {
        for (int j=0; j<ncols; j++)
        {
            p1 = nodes[connectivity(i,j)]->getPnt();
            if (j < ncols-1)
            {
                p2 = nodes[connectivity(i,j+1)]->getPnt();
            }
            else
            {
                p2 = nodes[connectivity(i,0)]->getPnt();
            }
            l = (p2-p1).norm();
            if (l < shortest)
            {
                shortest = l;
            }
        }
    }
    return shortest;
}

bool geometry::isLiftingSurf(int currentID, std::vector<int> wakeIDs)
{
    for (size_t i=0; i<wakeIDs.size(); i++)
    {
        if (wakeIDs[i]-10000 == currentID)
        {
            return true;
        }
    }
    return false;
}

void geometry::createSurfaces(const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> &connectivity, const Eigen::MatrixXd &norms, const Eigen::VectorXi &allID )
{
    surface* s = nullptr;
    wake* w = nullptr;
    bodyPanel* bPan;
    wakePanel* wPan;
    std::vector<edge*> pEdges;
    for (Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index i=0; i<static_cast<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index>(nTris); i++)
    {
        std::vector<cpNode*> pNodes;
        pNodes.push_back(nodes[connectivity(i,0)]);
        pNodes.push_back(nodes[connectivity(i,1)]);
        pNodes.push_back(nodes[connectivity(i,2)]);
        pEdges = panEdges(pNodes); //Create edge or find edge that already exists
        if (allID(i) <= 1000)
        {
            if (i==0 || allID(i) != allID(i-1))
            {
                s = new surface(allID(i),this);
                surfaces.push_back(s);
            }
            bPan = new bodyPanel(pNodes,pEdges,norms.row(i),s,static_cast<size_t>(allID(i)));
            s->addPanel(bPan);
        }
        else
        {
            if (i==0 || allID(i)!=allID(i-1))
            {
                w = new wake(static_cast<size_t>(allID(i)),this);
                wakes.push_back(w);
            }
            wPan = new wakePanel(pNodes,pEdges,norms.row(i),w,static_cast<size_t>(allID(i)));
            w->addPanel(wPan);
        }
    }
}

std::vector<edge*> geometry::panEdges(const std::vector<cpNode*>  &pNodes)
{
    size_t i1,i2;
    std::vector<edge*> triEdges;
    edge* e;
    for (size_t i=0; i<pNodes.size(); i++)
    {
        i1 = i;
        if (i == pNodes.size()-1)
        {
            i2 = 0;
        }
        else
        {
            i2 = i+1;
        }
        e = findEdge(pNodes[i1],pNodes[i2]);
        triEdges.push_back(e);
    }
    return triEdges;
}

edge* geometry::findEdge(cpNode* n1,cpNode* n2)
{
    for (edges_index_type i=0; i<edges.size(); i++)
    {
        if (edges[i]->sameEdge(n1, n2))
        {
            return edges[i];
        }
    }

    // If edge doesn't exist, create one
// NOTE: edge never used geometry that was linked to it
//    edge* e = new edge(n1,n2,this);
    edge* e = new edge(n1,n2);
    edges.push_back(e);
    return e;
}

void geometry::createOctree()
{
    std::vector<panel*> panels;
    std::vector<wakePanel*> tempW;
    std::vector<bodyPanel*> tempB;
    
    for (surfaces_index_type i=0; i<surfaces.size(); i++)
    {
        tempB = surfaces[i]->getPanels();
        panels.insert(panels.end(),tempB.begin(),tempB.end());
    }
    for (wakes_index_type i=0; i<wakes.size(); i++)
    {
        tempW = wakes[i]->getPanels();
        panels.insert(panels.end(),tempW.begin(),tempW.end());
    }
    pOctree.addData(panels);
}


void geometry::setInfCoeff()
{
	// Construct doublet and source influence coefficient matrices for body panels
	size_t nBodyPans = bPanels.size();
	size_t nWakePans = wPanels.size();
	size_t nPans = nBodyPans + nWakePans;

	Eigen::Matrix<size_t, 9, 1> percentage;
	percentage << 10, 20, 30, 40, 50, 60, 70, 80, 90;

	A.resize(static_cast<Eigen::MatrixXd::Index>(nBodyPans), static_cast<Eigen::MatrixXd::Index>(nBodyPans));
	B.resize(static_cast<Eigen::MatrixXd::Index>(nBodyPans), static_cast<Eigen::MatrixXd::Index>(nBodyPans));
	C.resize(static_cast<Eigen::MatrixXd::Index>(nBodyPans), static_cast<Eigen::MatrixXd::Index>(w2Panels.size()));

	A.setZero();
	B.setZero();
	C.setZero();

	for (size_t j = 0; j < nBodyPans; j++)
	{
		Eigen::Matrix<double, 1, Eigen::Dynamic> Arow = A.row(j);
		for (size_t i = 0; i < nBodyPans; i++)
		{
			bPanels[j]->panelPhiInf(bPanels[i]->getCenter(), B(static_cast<Eigen::MatrixXd::Index>(i), static_cast<Eigen::MatrixXd::Index>(j)), A(static_cast<Eigen::MatrixXd::Index>(i), static_cast<Eigen::MatrixXd::Index>(j)));
		}
		for (int i = 0; i < percentage.size(); i++)
		{
			if ((100 * j / nPans) <= percentage(i) && 100 * (j + 1) / nPans > percentage(i))
			{
				std::cout << percentage(i) << "%\t" << std::flush;
			}
		}
	}

	for (size_t i = 0; i<nBodyPans; i++)
	{
		bPanels[i]->setIndex(static_cast<int>(i));
	}

	std::vector<bodyPanel*> interpPans(4); // [Upper1 Lower1 Upper2 Lower2]  Panels that start the bounding wakelines of the wake panel.  Doublet strength is constant along wakelines (muUpper-muLower) and so the doublet strength used for influence of wake panel is interpolated between wakelines.
	double interpCoeff;
	double influence;

	Eigen::Vector4i indices;
	for (size_t j = 0; j < nWakePans; j++)
	{
		wPanels[j]->interpPanels(interpPans, interpCoeff);
		indices = interpIndices(interpPans);
		for (size_t i = 0; i < nBodyPans; i++)
		{
			Eigen::MatrixXd::Index ii(static_cast<Eigen::MatrixXd::Index>(i));
			influence = wPanels[j]->dubPhiInf(bPanels[i]->getCenter());
			A(ii, indices(0)) += influence * (1 - interpCoeff);
			A(ii, indices(1)) += influence * (interpCoeff - 1);
			A(ii, indices(2)) += influence * interpCoeff;
			A(ii, indices(3)) -= influence * interpCoeff;
		}
		for (int i = 0; i < percentage.size(); i++)
		{
			if ((100 * (nBodyPans + j) / nPans) <= percentage(i) && 100 * (nBodyPans + j + 1) / nPans > percentage(i))
			{
				std::cout << percentage(i) << "%\t" << std::flush;
			}
		}
	}

	// Construct doublet influence coefficient matrices for bufferWake panels

	for (size_t i = 0; i<nBodyPans; i++)
	{
		for (wakePanels_index_type j = 0; j<w2Panels.size(); j++)
		{
			C(static_cast<Eigen::MatrixXd::Index>(i), static_cast<Eigen::MatrixXd::Index>(j)) = w2Panels[j]->dubPhiInf(bPanels[i]->getCenter());
		}
	}

	std::cout << "Complete" << std::endl;

	if (writeCoeffFlag)
	{
		writeInfCoeff();
	}
}


void geometry::linSetInfCoeff()
{
	// Construct doublet and source influence coefficient matrices for body panels
	size_t nBodyPans = bPanels.size();
	size_t nWakePans = wPanels.size();
	size_t nPans = nBodyPans + nWakePans;

	size_t nBodyNodes = bodyNodes.size();

	Eigen::Matrix<size_t, 9, 1> percentage;
	percentage << 10, 20, 30, 40, 50, 60, 70, 80, 90;

	A.resize(static_cast<Eigen::MatrixXd::Index>(nBodyNodes), static_cast<Eigen::MatrixXd::Index>(nBodyNodes));
	B.resize(static_cast<Eigen::MatrixXd::Index>(nBodyNodes), static_cast<Eigen::MatrixXd::Index>(nBodyPans));
	A.setZero();
	B.setZero();

	Eigen::Matrix<double, 1, Eigen::Dynamic> Arow;
	Eigen::Vector3d ctrlPnt;

	for (size_t i = 0; i < nBodyPans; i++)
	{
		bPanels[i]->linTransformPanel();
	}

	for (size_t i = 0; i < nBodyNodes; i++)
	{
		ctrlPnt = nodes[i]->calcCP();
		Arow = A.row(i);
		for (size_t j = 0; j < nBodyPans; j++)
		{
			bPanels[j]->linDubPhiInf(ctrlPnt, Arow);
			bPanels[j]->srcPanelPhiInf(ctrlPnt, B(static_cast<Eigen::MatrixXd::Index>(i), static_cast<Eigen::MatrixXd::Index>(j)));
		}
		A.row(i) = Arow;

		for (int j = 0; j < percentage.size(); j++)
		{
			if ((100 * i / nBodyNodes) <= percentage(j) && 100 * (i + 1) / nBodyNodes > percentage(j))
			{
				std::cout << percentage(j) << "%\t" << std::flush;
			}
		}
	}

	for (size_t i = 0; i<nBodyPans; i++)
	{
		bPanels[i]->setIndex(static_cast<int>(i));
	}

	std::cout << "Complete" << std::endl;

	if (writeCoeffFlag)
	{
		writeInfCoeff();
	}
}


void geometry::supSetInfCoeff()
{
	// Construct doublet and source influence coefficient matrices for body panels

	// No wakes for the time being, so no 'subsonic' trailing edges allowed

	size_t nBodyPans = bPanels.size();
	size_t nWakePans = wPanels.size();
	size_t nPans = nBodyPans + nWakePans;

	size_t nBodyNodes = bodyNodes.size();

	Eigen::Matrix<size_t, 9, 1> percentage;
	percentage << 10, 20, 30, 40, 50, 60, 70, 80, 90;

	A.resize(static_cast<Eigen::MatrixXd::Index>(nBodyNodes), static_cast<Eigen::MatrixXd::Index>(nBodyNodes));
	B.resize(static_cast<Eigen::MatrixXd::Index>(nBodyNodes), static_cast<Eigen::MatrixXd::Index>(nBodyPans));
	A.setZero();
	B.setZero();

	double Bmach = sqrt(pow(inputMach, 2) - 1);
	Eigen::Matrix<double, 1, Eigen::Dynamic> Arow;
	Eigen::Vector3d ctrlPnt, windDir;
	bool DODflag; // Domain of Dependence flag: true if panel is inside or intersecting Mach cone of control point, false otherwise
	windDir = supComputeWindDir(); // Used in DoD checking

	// Compute transformation matrices, and perform the transformation, for each panel
	for (size_t i = 0; i < nBodyPans; i++)
	{
		// A lot of this only needs to be computed once. Could put outside loop and input into function to increase efficiency
		bPanels[i]->supTransformPanel(alpha, beta, inputMach);
	}

	// Build A and B matrices
	for (size_t i = 0; i < nBodyNodes; i++)
	{
		// Control point in original Ref CSYS. Transformed to panel local CSYS in supPhiInf()
		ctrlPnt = nodes[i]->calcCP();

		Arow = A.row(i);
		for (size_t j = 0; j < nBodyPans; j++)
		{
			DODflag = bPanels[j]->supDODcheck(ctrlPnt, inputMach, windDir);
			bPanels[j]->supPhiInf(ctrlPnt, Arow, B(static_cast<Eigen::MatrixXd::Index>(i), static_cast<Eigen::MatrixXd::Index>(j)), DODflag, inputMach);
			//bPanels[j]->supPhiInf11012019(ctrlPnt, Arow, B(static_cast<Eigen::MatrixXd::Index>(i), static_cast<Eigen::MatrixXd::Index>(j)), DODflag, inputMach, i);
		}
		A.row(i) = Arow;

		// Completion percentage (not exactly accurate w/ ss scheme, still useful though)
		for (int j = 0; j < percentage.size(); j++)
		{
			if ((100 * i / nBodyNodes) <= percentage(j) && 100 * (i + 1) / nBodyNodes > percentage(j))
			{
				std::cout << percentage(j) << "%\t" << std::flush;
			}
		}
	}

	for (size_t i = 0; i<nBodyPans; i++)
	{
		bPanels[i]->setIndex(static_cast<int>(i));
	}

	std::cout << "Complete" << std::endl;

	if (writeCoeffFlag)
	{
		writeInfCoeff();
	}
}


Eigen::Vector3d geometry::supComputeWindDir()
{
	Eigen::Vector3d windDir;
	double aalpha, bbeta;

	if (alpha == 0 && beta == 0)
	{
		windDir << 1, 0, 0;
	}
	else
	{
		aalpha = alpha * M_PI / 180;
		bbeta = beta * M_PI / 180;

		windDir << cos(aalpha)*cos(bbeta), -sin(bbeta), sin(aalpha)*cos(bbeta);
	}

	return windDir;
}


void geometry::setBodyToWind(double a, double b)
{
	a *= M_PI / 180;
	b *= M_PI / 180;

	body2wind << cos(a)*cos(b), -sin(b), sin(a)*cos(b),
		cos(a)*sin(b), cos(b), sin(a)*sin(b),
		-sin(a), 0, cos(a);
}


//Eigen::Vector3d geometry::bodyToWind(double vecIn, double a, double b)
//{
//	Eigen::Vector3d vecOut;
//
//	
//
//	vecOut = body2wind * vecIn;
//
//	return vecOut;
//}


Eigen::Vector4i geometry::interpIndices(std::vector<bodyPanel*> interpPans)
{
    Eigen::Vector4i indices;
    for (size_t i=0; i<interpPans.size(); i++)
    {
        indices(static_cast<Eigen::Vector4i::Index>(i)) = interpPans[i]->getIndex();
    }
    return indices;
}


std::vector<surface*> geometry::getSurfaces()
{
    return surfaces;
}

std::vector<wake*> geometry::getWakes()
{
    return wakes;
}

std::vector<panel*> geometry::getPanels()
{
    std::vector<panel*> panels;
    
    for (surfaces_index_type i=0; i<surfaces.size(); i++)
    {
        std::vector<bodyPanel*> temp = surfaces[i]->getPanels();
        panels.insert(panels.end(),temp.begin(),temp.end());
    }
    for (wakes_index_type i=0; i<wakes.size(); i++)
    {
        std::vector<wakePanel*> temp = wakes[i]->getPanels();
        panels.insert(panels.end(),temp.begin(),temp.end());
    }
    return panels;
}

bool geometry::infCoeffFileExists()
{
    boost::filesystem::path p = boost::filesystem::current_path().string()+"/" + infCoeffFile;
    if (boost::filesystem::exists(p))
    {
        std::ifstream fid;
        fid.open(p.string());
        size_t nPans;
        fid >> nPans;
        fid.close();
        if (nPans != bPanels.size())
        {
            return false;
        }
        return true;
    }

    return false;
}

void geometry::readInfCoeff()
{
    std::cout << "Reading Influence Coefficients from " << infCoeffFile << "..." << std::endl;

    std::ifstream fid;
    fid.open(infCoeffFile);
    int nPans, nW2pans;
    fid >> nPans;
    A.resize(nPans,nPans);
    B.resize(nPans,nPans);
    for (Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index i=0; i<static_cast<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index>(bPanels.size()); i++)
    {
        for (Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index j=0; j<static_cast<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index>(bPanels.size()); j++)
        {
            fid >> A(i,j);
        }
    }
    for (Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index i=0; i<static_cast<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index>(bPanels.size()); i++)
    {
        for (Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index j=0; j<static_cast<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index>(bPanels.size()); j++)
        {
            fid >> B(i,j);
        }
    }
    if (vortPartFlag)
    {
        fid >> nW2pans;
        C.resize(nPans, nW2pans);
        
        for (Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index i=0; i<static_cast<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index>(bPanels.size()); i++)
        {
            for (Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index j=0; j<static_cast<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index>(bPanels.size()); j++)
            {
                fid >> C(i,j);
            }
        }
        
    }
    fid.close();
}

void geometry::writeInfCoeff()
{
    std::cout << "Writing Influence Coefficients to " << infCoeffFile << "..." << std::endl;
    std::ofstream fid;
    fid.open(infCoeffFile);
    fid << bPanels.size() << "\n";
    for (Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index i=0; i<static_cast<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index>(bPanels.size()); i++)
    {
        for (Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index j=0; j<static_cast<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index>(bPanels.size()); j++)
        {
            fid << A(i,j) << "\t";
        }
        fid << "\n";
    }
    for (Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index i=0; i<static_cast<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index>(bPanels.size()); i++)
    {
        for (Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index j=0; j<static_cast<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index>(bPanels.size()); j++)
        {
            fid << B(i,j) << "\t";
        }
        fid << "\n";
    }
    if(vortPartFlag)
    {
        fid << w2Panels.size() << "\n";
        for (Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index i=0; i<static_cast<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index>(bPanels.size()); i++)
        {
            for (Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index j=0; j<static_cast<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>::Index>(bPanels.size()); j++)
            {
                fid << C(i,j) << "\t";
            }
            fid << "\n";
        }
    }
    fid.close();
}

Eigen::MatrixXd geometry::getNodePnts()
{
    Eigen::MatrixXd nodePnts(nodes.size(),3);
    for (nodes_index_type i=0; i<nodes.size(); i++)
    {
        nodePnts.row(static_cast<Eigen::MatrixXd::Index>(i)) = nodes[i]->getPnt();
    }
    return nodePnts;
}

double geometry::pntPotential(const Eigen::Vector3d &pnt, const Eigen::Vector3d &Vinf)
{
    double pot = 0;
    for (bodyPanels_index_type i=0; i<bPanels.size(); i++)
    {
        pot += bPanels[i]->panelPhi(pnt);
    }
    for (wakePanels_index_type i=0; i<wPanels.size(); i++)
    {
        pot += wPanels[i]->panelPhi(pnt);
    }
    pot += Vinf.dot(pnt);
    return pot;
}

double geometry::wakePotential(const Eigen::Vector3d &pnt)
{
    double pot = 0;
    for (wakePanels_index_type i=0; i<wPanels.size(); i++)
    {
        pot += wPanels[i]->panelPhi(pnt);
    }
    return pot;
}

Eigen::Vector3d geometry::pntVelocity(const Eigen::Vector3d &pnt, const Eigen::Vector3d &Vinf, double PG)
{
    Eigen::Vector3d vel = Eigen::Vector3d::Zero();
    for (bodyPanels_index_type i=0; i<bPanels.size(); i++)
    {
        vel += bPanels[i]->panelV(pnt);
    }
    for (wakePanels_index_type i=0; i<wPanels.size(); i++)
    {
        vel += wPanels[i]->panelV(pnt);
    }
    vel += Vinf;

    vel(0) /= PG; // Prandtl-Glauert Correction
    return vel;
}


void geometry::clusterCheck()
{
    std::ofstream fid;
    fid.open("ClusterCheck.txt");
    int index;
    std::vector<bodyPanel*> clust;
    for (bodyPanels_index_type i=0; i<bPanels.size(); i++)
    {
        clust = bPanels[i]->getCluster();
        for (size_t j=0; j<clust.size(); j++)
        {
            index = (int)std::distance(bPanels.begin(),std::find(bPanels.begin(), bPanels.end(), clust[j]));
            fid << index+1 << "\t";
        }
        fid << "\n";
    }
    fid.close();
}

void geometry::calcTimeStep(){
    
    // Timestep will be set so that the step in the streamwise direction results in the same distance as the particles are spaced apart to allow for equal particle spacing and thus sufficient overlap. Sized for the average wake panel width becuase panel clustering at geometry tips doesn't represent average spacing.
    
    
    if(dt == 0) // If timestep is set in input file, it will not be zero and won't be modified
    {
        std::vector<edge*> Tedges;
        for(edges_index_type i=0; i<edges.size(); i++)
        {
            if(edges[i]->isTE())
            {
                Tedges.push_back(edges[i]);
            }
        }
        
        double sum = 0;
        for(size_t i=0; i < Tedges.size(); i++)
        {
            sum += Tedges[i]->length();
        }
        double avgLength = sum/static_cast<double>(Tedges.size());
        
        dt = avgLength/inputV;
        
        std::cout << "\tCalculated time step : " << dt << " sec" << std::endl;
    }
    
}


void geometry::createVPWakeSurfaces(const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> &wakeConnectivity, const Eigen::MatrixXd &wakeNorms,  const std::vector<size_t> &VPwakeID,  std::vector<bool> iisFirstPanel){
    
    wake* w = nullptr;
    wakePanel* wPan;
    wakePanel* firstPan; // Used only to keep track of panel for setting bw2's parent
    std::vector<edge*> pEdges;
    
    w = new wake(VPwakeID[0],this);
    wakes.push_back(w);
    
    for (int i=0; i<wakeConnectivity.rows(); i++)
    {
        std::vector<cpNode*> pNodes;
        for(int j=0; j<wakeConnectivity.row(i).size(); j++)
        {
            pNodes.push_back(nodes[wakeConnectivity(i,j)]);
        }
        pEdges = panEdges(pNodes); // Create edges or find edges that already exists
        
        
        wPan = new wakePanel(pNodes,pEdges,wakeNorms.row(i),w,VPwakeID[static_cast<size_t>(i)]);
        w->addPanel(wPan);
        
        if( iisFirstPanel[static_cast<size_t>(i)] )
        {
            firstPan = wPan;
            wPan->isSecondRow = false;
        }
        else
        {
            wPan->setBufferParent(firstPan); // There will ALWAYS be a 'firstPanel' before second row
            wPan->isSecondRow = true;
            w2Panels.push_back(wPan);
        }
    }
}




//void geometry::moveGeom( std::vector<double> bodyKin ){
//    
//    
//    for (int i=0; i<nodes.size(); i++) {
//        
//        Eigen::Vector3d localMovement;
//        
//        Eigen::Vector3d pos = nodes[i]->getPnt();
//        
//        // U = U3 + (-q*z + r*y)
//        localMovement.x() = bodyKin[0] - bodyKin[4]*pos.z() + bodyKin[5]*pos.y();
//        
//        // V = V3 + (-r*x + p*z)
//        localMovement.y() = bodyKin[1] - bodyKin[5]*pos.x() + bodyKin[3]*pos.z();
//        
//        // W = W3 + (-p*y + q*x)
//        localMovement.z() = bodyKin[2] - bodyKin[3]*pos.y() + bodyKin[4]*pos.x();
//        
//        nodes[i]->setPnt( pos + localMovement*dt );
//    }
//    
//    // Re-calculate the panel geometry (center, normal, etc.)
//    for (int i=0; i<bPanels.size(); i++) {
//        bPanels[i]->setGeom();
//    }
//    
//    for(int i=0; i<wPanels.size(); i++){
//        wPanels[i]->setGeom();
//    }
//    
//    for(int i=0; i<w2Panels.size(); i++){
//        w2Panels[i]->setGeom();
//    }
//    
//    std::cout << "Geom moved..." << std::endl;
//    
//}



//bool linDubFlag = false;
//if (inputMach > 1)
//{
//	linDubFlag = true;
//}

//for (size_t j = 0; j<nBodyPans; j++)
//{
//	for (size_t i = 0; i<nDubBodyCPs; i++)
//	{
//		bPanels[j]->panelPhiInf(bDubCPs[i], A(static_cast<Eigen::MatrixXd::Index>(i), static_cast<Eigen::MatrixXd::Index>(j)), true, linDubFlag);
//	}
//	for (size_t i = 0; i<nSrcBodyCPs; i++)
//	{
//		bPanels[j]->panelPhiInf(bSrcCPs[i], B(static_cast<Eigen::MatrixXd::Index>(i), static_cast<Eigen::MatrixXd::Index>(j)), false, linDubFlag);
//	}
//	/////////////////////////////////// UPDATE //////////////////////////////////////
//	for (int i = 0; i<percentage.size(); i++)
//	{
//		if ((100 * j / nPans) <= percentage(i) && 100 * (j + 1) / nPans > percentage(i))
//		{
//			std::cout << percentage(i) << "%\t" << std::flush;
//		}
//	}
//}


//Eigen::Vector3d nodePnt;
//Eigen::Vector3d tempNorm;
//tempNorm.setZero();

//// WRITE TO
//std::ofstream fid;
//std::string myFile = "cntrlPnts.csv";
//fid.open(myFile);
//fid << "nodes_x" << "," << "nodes_y" << "," << "nodes_z" << ",";
//fid << "norm_x" << "," << "norm_y" << "," << "norm_z" << ",";
//fid << "CP_x" << "," << "CP_y" << "," << "CP_z";
//fid << "\n";

//if (inMach > 1)
//{
//	// WRITE TO
//	size_t k = 0;
//	size_t m = 0;

//	for (nodes_index_type i = 0; i < nodes.size(); i++)
//	{
//		tempNorm.setZero();
//		nodePnt = nodes[i]->getPnt();

//		// check if node is part of a body panel
//		if (nodes[i]->getBodyPans().size() > 0)
//		{
//			for (bodyPanels_index_type j = 0; j < nodes[i]->getBodyPans().size(); j++)
//			{
//				tempNorm += nodes[i]->getBodyPans()[j]->getNormal();
//			}
//			tempNorm.normalize();
//			bodyCPs.push_back(nodes[i]->getPnt() - scaleNorm * tempNorm);

//			// WRITE TO
//			fid << nodePnt(0) << "," << nodePnt(1) << "," << nodePnt(2) << ",";
//			fid << tempNorm(0) << "," << tempNorm(1) << "," << tempNorm(2) << ",";
//			fid << bodyCPs[k](0) << "," << bodyCPs[k](1) << "," << bodyCPs[k](2);
//			++k;
//		}
//		// if node is part of wake panel, CP is the node
//		else
//		{
//			wakeCPs.push_back(nodes[i]->getPnt());

//			// WRITE TO
//			fid << wakeCPs[m](0) << "," << wakeCPs[m](1) << "," << wakeCPs[m](2) << ",";
//			fid << wakeCPs[m](0) << "," << wakeCPs[m](1) << "," << wakeCPs[m](2);
//			++m;
//		}
//		// WRITE TO
//		fid << "\n";
//	}
//}
//else
//{
//	for (bodyPanels_index_type i = 0; i < bPanels.size(); i++)
//	{
//		bodyCPs.push_back(bPanels[i]->getCenter());
//	}
//	for (wakePanels_index_type i = 0; i < wPanels.size(); i++)
//	{
//		wakeCPs.push_back(wPanels[i]->getCenter());
//	}
//}
//// WRITE TO
//fid.close();