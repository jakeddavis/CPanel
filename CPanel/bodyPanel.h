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
 *    David D. Marshall - misc. changes
 ******************************************************************************/

#ifndef __CPanel__bodyPanel__
#define __CPanel__bodyPanel__

#include <iostream>
#include <math.h>
#include "panel.h"
#include "particle.h"
//#include "edge.h"
//#include "cpNode.h"
//
//class edge;
class cpNode;
class particle;

class bodyPanel : public panel
{
	using bodyPanels_type = std::vector<bodyPanel *>;
	using bodyPanels_index_type = bodyPanels_type::size_type;

	/*using nodes_type = std::vector<cpNode *>;
	using nodes_index_type = nodes_type::size_type;*/

    surface* parentSurf;
    bodyPanels_type neighbors;
    bodyPanels_type cluster;
    int TSorder;
    double sourceStrength = 0; // Initialize with zero so can print before compVelocity
    bool upper; // Sheds wake panel from lower edge
    bool lower; // Sheds wake panel from upper edge
    bool TEpanel;
    edge* TE;
    bool tipFlag;
    bool streamFlag; // Surface Streamline crosses panel.
    Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
    double Cp = 0;
    
    int index; // Index in panel vector contained in geometry class.  Used for interpolating strength for wake panel influences.
    
    double srcSidePhi(const double &PN,const double &Al, const double &phiV,const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s);
    Eigen::Vector3d srcSideV(const double &PN,const double &Al,const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s,const Eigen::Vector3d &l,const Eigen::Vector3d &m,const Eigen::Vector3d &n);
    inline double pntSrcPhi(const double &PJK);
    inline Eigen::Vector3d pntSrcV(const Eigen::Vector3d &pjk);
    
    Eigen::Vector3d velocity2D(const Eigen::Vector3d &pnt,double pntPotential);
    Eigen::Vector3d velocity3D(const Eigen::Vector3d &pnt,double pntPotential);
    
    bool clusterTest(bodyPanel* other, double angle,bool upFlag,bool lowFlag);
    bool nearTrailingEdge();
    
    Eigen::Matrix3d velocityGradientPointSource(Eigen::Vector3d POI);
    Eigen::Matrix3d velocityGradientQuadSource(Eigen::Vector3d POI);
    Eigen::Matrix3d velocityGradientTriSource(Eigen::Vector3d POI);

	//std::vector<bool> edgeFlags;

	Eigen::Matrix3d supTransMat;
	std::vector<Eigen::Vector3d> supLocalNodes;
	double supAreaCorrect;
    
public:
    bodyPanel(std::vector<cpNode*> nodes, std::vector<edge*> pEdges, Eigen::Vector3d bezNorm,surface* parentSurf, size_t surfID);
        
    void addNeighbor(bodyPanel* p);
    void setUpper();
    void setLower();
    void setTEpanel(edge* trailingEdge);
    void setIndex(int i);
    void setLSflag();
    void setTipFlag();
    void setCluster();
    
    double panelPhi(const Eigen::Vector3d &POI);
    Eigen::Vector3d panelV(const Eigen::Vector3d &POI);

    void panelPhiInf(const Eigen::Vector3d &POI, double &phiSrc,double &phiDub);
    void panelVInf(const Eigen::Vector3d &POI, Eigen::Vector3d &vSrc,Eigen::Vector3d &vDub);
    Eigen::Vector3d pntVInf(const Eigen::Vector3d &POI);
    
    Eigen::Vector3d pntVelocity(const Eigen::Vector3d &pnt,double pntPotential, double PG, const Eigen::Vector3d &Vinf);
    
    double dist2Pan(bodyPanel* other);
    
    void computeVelocity(double PG, const Eigen::Vector3d &Vinf);
    void computeCp(double Vinf);
    void computeCp(double Vinf, double dt);
    
    Eigen::Vector3d computeMoments(const Eigen::Vector3d &cg);
    
    void setSigma(Eigen::Vector3d Vinf, double Vnorm);
    void setMu(double dubStrength);
    
    void setStreamFlag();
    
    std::vector<bodyPanel*> getNeighbors() {return neighbors;}
    std::vector<bodyPanel*> getRelatedPanels();
    
    double getSigma() {return sourceStrength;}
    double getMu() {return doubletStrength;}
    bool isUpper() {return upper;}
    bool isLower() {return lower;}
    bool isTEpanel() {return TEpanel;}
    edge* getTrailingEdge() {return TE;}
    
    bool isTipPan() {return tipFlag;}
    bool getStreamFlag() {return streamFlag;}
    int getIndex() {return index;}
    Eigen::Vector3d getGlobalV() {return velocity;}
    double getCp() {return Cp;}
    
    std::vector<bodyPanel*> getCluster() {return cluster;}
    
    Eigen::Vector3d partStretching(particle* part);

	// ss
	void srcPanelPhiInf(const Eigen::Vector3d &POI, double &phi);
	void linComputeVelocity(double PG, Eigen::Vector3d &Vinf);
	Eigen::Vector3d linComputeVelocity2(double PG, Eigen::Vector3d &Vinf, Eigen::Vector3d &POI);
	double linGetTEdubStrength();

	void supPhiInf(const Eigen::Vector3d &P, Eigen::Matrix<double, 1, Eigen::Dynamic> &Arow, double &Phi, bool DOIflag, const double mach, Eigen::Vector3d &windDir);
	bool supDOIcheck(Eigen::Vector3d &POI,const double Mach, Eigen::Vector3d &windDir);
	//bool supEdgeCheck(edge* myEdge, Eigen::Vector3d &P, const double B);
	void supTransformPanel(const double Bmach, double alpha, double beta, const double M);
	//Eigen::Matrix3d supG2LSmatrix(const double Bmach);
	void supSetG2LSmatrix(const double Bmach, const double a, const double b, const double M);
	Eigen::Vector3d supConePanelInter(const Eigen::Vector3d &POI, const double Mach, Eigen::Vector3d &windDir);
	Eigen::Matrix3d supGetLocalSys(Eigen::Vector3d &windDir);

	Eigen::Vector2d supEdgeInfSon(const double ym1, const double ym2, const double xmc, const double ym1c, const double ym2c, const double R1, const double R2, const double lam, const double z);
	Eigen::Vector2d supEdgeInfSub(const double R1, const double R2, const double ym1c, const double ym2c, const double xmc, const double m, const double z, const double eps1, const double eps2, bool mFlag);
	Eigen::Vector2d supEdgeInfSup(const double R1, const double R2, const double ym1, const double ym2, const double xm, const double lam, const double z, const double eps1, const double eps2);

	void supOutputGeom(const Eigen::Vector3d &POI, bool outPOI);

	void supComputeVelocity(Eigen::Vector3d &Vinf);
	void supComputeCp();

	//std::vector<bool> getEdgeFlags() { return edgeFlags; }
    
};

#endif /* defined(__CPanel__bodyPanel__) */
