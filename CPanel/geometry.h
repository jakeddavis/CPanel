/*******************************************************************************
 * Copyright (c) 2015 Chris Satterwhite
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

//#ifndef __CPanel__geometry__
//#define __CPanel__geometry__
//
//#include <iostream>
//#include <vector>
//#include <Eigen/Dense>
//#include <fstream>
//#include <array>
//#include <boost/filesystem/operations.hpp>
//#include <boost/filesystem/path.hpp>
//#include "panelOctree.h"
//#include "surface.h"
//#include "liftingSurf.h"
//#include "wake.h"
//#include "bodyPanel.h"
//#include "wakePanel.h"
//#include "edge.h"
//#include "inputParams.h"
//#include "cpNode.h"
////#include "octreeFile.h" //Taking this out just for CPanel test. Need to put back in after I figure it out...
//
//class geometry
//{
//    std::vector<surface*> surfaces;
//    std::vector<wake*> wakes;
////    std::vector<liftingSurf*> liftingSurfs;
////    std::vector<surface*> nonLiftingSurfs;
//    std::vector<bodyPanel*> bPanels;
//    std::vector<wakePanel*> wPanels;
//    std::vector<wakePanel*> w2Panels; // Buffer wake row two
//    double c_w;
//    double inputV;
//
//    panelOctree pOctree;
//    std::vector<cpNode*> nodes;
//    std::vector<edge*> edges;
//    short nNodes;
//    short nTris;
//    double dt; //making public so cpCase can access it
//
//
//    Eigen::MatrixXd B; // Source Influence Coefficient Matrix
//    Eigen::MatrixXd A; // Doublet Influence Coefficient Matrix
//
//    Eigen::MatrixXd C; // Wake Doublet Influence Coefficient Matrix
//
//    bool writeCoeffFlag;
//    std::string infCoeffFile;
//    bool vortPartFlag;
//    std::vector<bool> isFirstPanel;
//
//    void readTri(std::string tri_file, bool normFlag);
//    std::vector<edge*> panEdges(const std::vector<cpNode*> &pNodes);
//    edge* findEdge(cpNode* n1,cpNode* n2);
//    void createSurfaces(const Eigen::MatrixXi &connectivity, const Eigen::MatrixXd &norms, const Eigen::VectorXi &allID, std::vector<int> wakeIDs, bool VortPartFlag);
//    void createVPWakeSurfaces(const Eigen::MatrixXi &wakeConnectivity, const Eigen::MatrixXd &wakeNorms,  const std::vector<int> &VPwakeID, std::vector<bool> isFirstPanel);
//
//    void createOctree();
////    void setTEPanels();
////    void setTEnodes();
//    void getLiftingSurfs(std::vector<surface*>& wakes, std::vector<surface*>& liftingSurfs);
//    void setNeighbors(panel* p,int targetID);
////    void scanNode(panel* p, node<panel>* current, node<panel>* exception);
//    bool isLiftingSurf(int currentID, std::vector<int> wakeIDs);
////    void correctWakeNodes(int wakeNodeStart);
//    void correctWakeConnectivity(int wakeNodeStart,int wakeTriStart,Eigen::MatrixXi &connectivity);
//    double shortestEdge(const Eigen::MatrixXi &connectivity);
//    liftingSurf* getParentSurf(int wakeID);
//
//    void setInfCoeff();
//
//    Eigen::Vector4i interpIndices(std::vector<bodyPanel*> interpPans);
//
//    bool infCoeffFileExists();
//    void readInfCoeff();
//    void writeInfCoeff();
//
//    bool unsteadySim;
//
//    std::string bodyKinFile;
//    Eigen::VectorXd bodyKin;
//
//    void readBodyKinFile();
//    Eigen::Vector3d Vinf(Eigen::Vector3d POI);
//
//
//public:
//
//    geometry(inputParams* p)
//    {
//        std::stringstream temp;
//        temp << p->geomFile->name << ".infCoeff";
//        infCoeffFile = temp.str();
//        writeCoeffFlag = p->writeCoeffFlag;
//        vortPartFlag = p->vortexParticles;
//        inputV = p->velocities(0);
//        dt = p->timeStep;
//        unsteadySim = p->unsteady;
//        if(unsteadySim){
//            bodyKinFile = p->bodyKinFileLoc;
//            readBodyKinFile();
//        }
//        readTri(p->geomFile->file, p->normFlag);
//    }
//
//    virtual ~geometry();
//
//    geometry(const geometry& copy);
//
//    geometry& operator=(const geometry &rhs);
//
//    double getDt() { return dt;}
//    double pntPotential(const Eigen::Vector3d &pnt, const Eigen::Vector3d &Vinf);
//    double wakePotential(const Eigen::Vector3d &pnt);
//    Eigen::Vector3d pntVelocity(const Eigen::Vector3d &pnt, const Eigen::Vector3d &Vinf, double PG);
//
//    void clusterCheck();
//
//    short getNumberOfNodes() {return nNodes;}
//    short getNumberOfTris() {return nTris;}
//    std::vector<cpNode*> getNodes() {return nodes;}
//    Eigen::MatrixXd getNodePnts();
////    std::vector<liftingSurf*> getLiftingSurfs() {return liftingSurfs;}
////    std::vector<surface*> getNonLiftingSurfs() {return nonLiftingSurfs;}
//
//    std::vector<surface*> getSurfaces();
//    panelOctree* getOctree() {return &pOctree;}
//    std::vector<panel*> getPanels();
//    std::vector<bodyPanel*>* getBodyPanels() {return &bPanels;}
//    std::vector<wakePanel*>* getWakePanels() {return &wPanels;}
//    std::vector<wakePanel*>* getBufferWake2Panels() {return &w2Panels;} // bw2
//    std::vector<wake*> getWakes();
//    Eigen::MatrixXd* getA() {return &A;}
//    Eigen::MatrixXd* getB() {return &B;}
//    Eigen::MatrixXd* getC() {return &C;} // bw2
//
//
//};
//
//
//#endif /* defined(__CPanel__geometry__) */



//
//  geometry.h
//  CPanel
//
//  Created by Chris Satterwhite on 4/30/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__geometry__
#define __CPanel__geometry__

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include <array>
#include "panelOctree.h"
#include "surface.h"
#include "liftingSurf.h"
#include "wake.h"
#include "bodyPanel.h"
#include "wakePanel.h"
#include "edge.h"
#include "inputParams.h"
#include "cpNode.h"

//#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

class geometry
{
	using surfaces_type = std::vector<surface *>;
	using surfaces_index_type = surfaces_type::size_type;
	using wakes_type = std::vector<wake *>;
	using wakes_index_type = wakes_type::size_type;
	using bodyPanels_type = std::vector<bodyPanel *>;
	using bodyPanels_index_type = bodyPanels_type::size_type;
	using wakePanels_type = std::vector<wakePanel *>;
	using wakePanels_index_type = wakePanels_type::size_type;

    surfaces_type surfaces;
    wakes_type wakes;
    bodyPanels_type bPanels;
    wakePanels_type wPanels;
    wakePanels_type w2Panels; // Buffer wake row two

    std::vector<bool> isFirstPanel;

    using nodes_type = std::vector<cpNode *>;
    using nodes_index_type = nodes_type::size_type;
    using edges_type = std::vector<edge *>;
    using edges_index_type = edges_type::size_type;

    panelOctree pOctree;
    nodes_type nodes;
    edges_type edges;
    size_t nNodes;
    size_t nTris;

	nodes_type bodyNodes;
	nodes_type wakeNodes;

    Eigen::MatrixXd A; // Doublet Influence Coefficient Matrix
    Eigen::MatrixXd B; // Source Influence Coefficient Matrix
    Eigen::MatrixXd C; // 2nd Row Buffer Wake Doublet Influence Coefficient Matrix



    bool writeCoeffFlag;
    bool vortPartFlag;
    std::string infCoeffFile;
    double dt;
    double inputV;

	// lin & ss
	bool subHOMFlag;
	double inputMach;
	double alpha;
	double beta;
	Eigen::Matrix3d body2wind;

    void readTri(std::string tri_file, bool normFlag);
    std::vector<edge*> panEdges(const std::vector<cpNode*> &pNodes);
    edge* findEdge(cpNode* n1,cpNode* n2);
    void createSurfaces(const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> &connectivity, const Eigen::MatrixXd &norms, const Eigen::VectorXi &allID );
    void createOctree();
    void getLiftingSurfs(std::vector<surface*>& wakes, std::vector<surface*>& liftingSurfs);
    void setNeighbors(panel* p,int targetID);
    bool isLiftingSurf(int currentID, std::vector<int> wakeIDs);
    void correctWakeConnectivity(size_t wakeNodeStart,size_t wakeTriStart,Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> &connectivity);
    double shortestEdge(const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> &connectivity);
    liftingSurf* getParentSurf(int wakeID);

    void setInfCoeff();
	void linSetInfCoeff();
	void supSetInfCoeff();
    Eigen::Vector4i interpIndices(std::vector<bodyPanel*> interpPans);

    bool infCoeffFileExists();
    void readInfCoeff();
    void writeInfCoeff();

    void createVPWakeSurfaces(const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> &wakeConnectivity, const Eigen::MatrixXd &wakeNorms,  const std::vector<size_t> &VPwakeID, std::vector<bool> isFirstPanel);
    void calcTimeStep();


	void supSortUpperLower(std::string tri_file);




public:
    geometry(inputParams* p)
    {
        std::stringstream temp;
        temp << p->geomFile->name << ".infCoeff";
        infCoeffFile = temp.str();
        writeCoeffFlag = p->writeCoeffFlag;
        vortPartFlag = p->vortexParticles;
        dt = p->timeStep;
        inputV = p->velocities(0);
        nNodes=0;
        nTris=0;

		subHOMFlag = p->subHOMFlag;
		inputMach = p->machs(0);
		alpha = p->alphas(0);
		beta = p->betas(0);


		//---------------------------------------------------------------------------------------//
		//---------------------------------------------------------------------------------------//

		// Read in supersonic sharp leading edge flag given by user (defined here for now)
		bool supSharpLEflag = true;
		// Check flag
		if (supSharpLEflag)
		{
			// Create new .tri file to write sorted data to
			std::string sorted = "/sorted.tri";
			std::string newSorted = p->geomFile->path + sorted;

			if (boost::filesystem::exists(newSorted))
			{
				boost::filesystem::remove(newSorted);
			}
				

			boost::filesystem::copy_file(p->geomFile->file, newSorted);

			// Sort tri file into upper and lower, and make copies of LE/TE nodes and edges
			supSortUpperLower(newSorted);


			/*std::cout << p->geomFile->path << std::endl;
			std::cout << p->geomFile->file << std::endl;*/
		}

		//---------------------------------------------------------------------------------------//
		//---------------------------------------------------------------------------------------//


        readTri(p->geomFile->file, p->normFlag);
    }

    virtual ~geometry();

    geometry(const geometry& copy);

    geometry& operator=(const geometry &rhs);

    double pntPotential(const Eigen::Vector3d &pnt, const Eigen::Vector3d &Vinf);
    double wakePotential(const Eigen::Vector3d &pnt);
    Eigen::Vector3d pntVelocity(const Eigen::Vector3d &pnt, const Eigen::Vector3d &Vinf, double PG);

    void clusterCheck();

    size_t getNumberOfNodes() {return nNodes;}
    size_t getNumberOfTris() {return nTris;}
    std::vector<cpNode*> getNodes() {return nodes;}
    Eigen::MatrixXd getNodePnts();

    std::vector<surface*> getSurfaces();
    panelOctree* getOctree() {return &pOctree;}
    std::vector<panel*> getPanels();
    std::vector<bodyPanel*>* getBodyPanels() {return &bPanels;}
    std::vector<wakePanel*>* getWakePanels() {return &wPanels;}
    std::vector<wakePanel*>* getBufferWake2Panels() {return &w2Panels;}
    std::vector<wake*> getWakes();
    Eigen::MatrixXd* getA() {return &A;}
    Eigen::MatrixXd* getB() {return &B;}
    Eigen::MatrixXd* getC() {return &C;}
    double getDt() {return dt;}

    void moveGeom( std::vector<double> bodyKin );

	std::vector<cpNode*> getBodyNodes() { return bodyNodes; }
	void geometry::setBodyToWind(double a, double b);
	Eigen::Vector3d geometry::supComputeWindDir();
};


#endif /* defined(__CPanel__geometry__) */
