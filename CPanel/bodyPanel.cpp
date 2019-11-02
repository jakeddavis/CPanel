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

#include "bodyPanel.h"
#include "edge.h"
#include "cpNode.h"
#include "surface.h"

#define _USE_MATH_DEFINES
#include <math.h>

bodyPanel::bodyPanel(std::vector<cpNode*> nnodes, std::vector<edge*> ppEdges,
		            Eigen::Vector3d bezNorm,surface* pparentSurf, size_t surfID)
  : panel(nnodes,ppEdges,bezNorm,surfID), parentSurf(pparentSurf), TSorder(3), upper(false), lower(false),
	TEpanel(false), TE(nullptr), tipFlag(false), streamFlag(false), index(-1)
{
    for (size_t i=0; i<pEdges.size(); i++)
    {
        pEdges[i]->addBodyPan(this);
    }
    for (size_t i=0; i<nodes.size(); i++)
    {
        nodes[i]->addBodyPanel(this);
    }
}

void bodyPanel::addNeighbor(bodyPanel* p)
{
    neighbors.push_back(p);
}

void bodyPanel::setUpper() {upper = true;}
void bodyPanel::setLower() {lower = true;}

void bodyPanel::setTEpanel(edge* trailingEdge)
{
    TEpanel = true;
    parentSurf->setLSflag();
    TE = trailingEdge;
}

void bodyPanel::setIndex(int i) {index = i;}

void bodyPanel::setTipFlag()
{
    if (!parentSurf->isLiftingSurf())
    {
        tipFlag = false;
        return;
    }
    
    int count = 0;
    double angle;
    Eigen::Vector3d nNormal;
    for (bodyPanels_index_type i=0; i<neighbors.size(); i++)
    {
        nNormal = neighbors[i]->getNormal();
        double dot = normal.dot(nNormal)/(normal.norm()*nNormal.norm());
        if (dot > 1)
        {
            dot = 1;
        }
        else if (dot < -1)
        {
            dot = -1;
        }
        
        angle = acos(dot);
        
        if (angle < pow(10,-8) && asin(nNormal(0)) > -M_PI/12)
        {
            // asin(nNormal(0)) is the angle normal vector makes with xy plane.  Catches bug where some leading edge panels were being flagged as on the tip patch.
            count++;
        }
    }
    if (count >= 2)
    {
        tipFlag = true;
    }
    
}

void bodyPanel::setSigma(Eigen::Vector3d Vinf, double Vnorm)
{
    sourceStrength = (-Vinf.dot(normal)+Vnorm); // CS: Vnorm is eminating from the panel
}

void bodyPanel::setMu(double dubStrength)
{
    doubletStrength = dubStrength;
}

void bodyPanel::setStreamFlag()
{
    streamFlag = true;
}

double bodyPanel::panelPhi(const Eigen::Vector3d &POI)
{
    double phiSrc,phiDub;
    phiSrc = 0;
    phiDub = 0;

    panelPhiInf(POI,phiSrc,phiDub);
    return -sourceStrength*phiSrc-doubletStrength*phiDub;
}

Eigen::Vector3d bodyPanel::panelV(const Eigen::Vector3d &POI)
{
    Eigen::Vector3d vSrc = Eigen::Vector3d::Zero();
    Eigen::Vector3d vDub = Eigen::Vector3d::Zero();
    
    panelVInf(POI,vSrc,vDub);
    
    return sourceStrength*vSrc+doubletStrength*vDub;
}

void bodyPanel::panelPhiInf(const Eigen::Vector3d &POI, double &phiSrc,double &phiDub)
{
    Eigen::Vector3d pjk = POI-center;
    bool itselfFlag = false;
    if (pjk.norm() < pow(10,-10))
    {
        itselfFlag = true;
    }
    Eigen::Matrix3d local = getLocalSys();
    double PN = pjk.dot(local.row(2));
    if (pjk.norm()/longSide > 5)
    {
        phiSrc = pntSrcPhi(pjk.norm());
        phiDub = pntDubPhi(PN,pjk.norm());
    }
    else
    {
		double Al;
		double phiV = 0;
        Eigen::Vector3d a,b,s;
        for (nodes_index_type i=0; i<nodes.size(); i++)
        {
            Eigen::Vector3d p1;
            Eigen::Vector3d p2;
            if (i!=nodes.size()-1)
            {
                p1 = nodes[i]->getPnt();
                p2 = nodes[i+1]->getPnt();
            }
            else
            {
                p1 = nodes[i]->getPnt();
                p2 = nodes[0]->getPnt();
            }
            a = POI-p1;
            b = POI-p2;
            s = p2-p1;
            Al = local.row(2).dot(s.cross(a));
            if (!itselfFlag)
            {
// Note: last parameter was not used
//                phiV = vortexPhi(PN,Al,a,b,s,local.row(0),local.row(1),local.row(2));
                phiV = vortexPhi(PN,Al,a,b,s,local.row(0),local.row(1));
                phiDub += phiV;
            }
            phiSrc += srcSidePhi(PN,Al,phiV,a,b,s);
            
        }
        phiSrc /= (4*M_PI);
        if (!itselfFlag)
        {
            phiDub /= (4*M_PI);
        }
        else
        {
            phiDub = -0.5;
        }
    }
}

void bodyPanel::srcPanelPhiInf(const Eigen::Vector3d &POI, double &phi)
{
	double a, m, d, e1, e2, h1, h2, r1, r2, nuXi, nuEta;
	double Q, J;

	Eigen::Vector3d pjk, POIloc;
	pjk = POI - center;
	POIloc = global2local(POI, true);

	double x, y, z, x1, y1, x2, y2;
	x = POIloc.x();
	y = POIloc.y();
	z = POIloc.z();

	if (pjk.norm() / longSide > 5)
	{
		phi = pntSrcPhi(pjk.norm());
	}
	else
	{
		for (nodes_index_type i = 0; i<nodes.size(); i++)
		{
			Eigen::Vector3d p1;
			Eigen::Vector3d p2;
			if (i != nodes.size() - 1)
			{
				p1 = linLocalNodes[i];
				p2 = linLocalNodes[i + 1];
			}
			else
			{
				p1 = linLocalNodes[i];
				p2 = linLocalNodes[0];
			}

			x1 = p1.x();
			y1 = p1.y();
			x2 = p2.x();
			y2 = p2.y();

			d = (p2 - p1).norm();
			nuEta = (x2 - x1) / d;
			nuXi = (y2 - y1) / d;
			a = (x - x1)*nuXi - (y - y1)*nuEta;
			m = (y2 - y1) / (x2 - x1);
			e1 = pow(x - x1,2) + pow(z,2);
			e2 = pow(x - x2,2) + pow(z,2);
			h1 = (x - x1) * (y - y1);
			h2 = (x - x2) * (y - y2);
			r1 = (POIloc - p1).norm();
			r2 = (POIloc - p2).norm();

			Q = log((r1 + r2 + d) / (r1 + r2 - d));
			J = atan((m*e1 - h1) / (z*r1)) - atan((m*e2 - h2) / (z*r2));

			phi += a * Q - abs(z)*J;

		}
		phi /= (4 * M_PI);
	}
}

//void bodyPanel::srcPanelPhiInf(const Eigen::Vector3d &POI, double &phi)
//{
//	supOutputGeom(POI, true);
//
//	Eigen::Vector3d pjk = POI - center;
//	bool itselfFlag = false;
//	if (pjk.norm() < pow(10, -10))
//	{
//		itselfFlag = true;
//	}
//	Eigen::Matrix3d local = getLocalSys();
//	double PN = pjk.dot(local.row(2));
//	if (pjk.norm() / longSide > 5)
//	{
//		phi = pntSrcPhi(pjk.norm());
//	}
//	else
//	{
//		double Al;
//		double phiV = 0;
//		Eigen::Vector3d a, b, s;
//		for (nodes_index_type i = 0; i<nodes.size(); i++)
//		{
//			Eigen::Vector3d p1;
//			Eigen::Vector3d p2;
//			if (i != nodes.size() - 1)
//			{
//				p1 = nodes[i]->getPnt();
//				p2 = nodes[i + 1]->getPnt();
//			}
//			else
//			{
//				p1 = nodes[i]->getPnt();
//				p2 = nodes[0]->getPnt();
//			}
//			a = POI - p1;
//			b = POI - p2;
//			s = p2 - p1;
//			Al = local.row(2).dot(s.cross(a));
//			if (!itselfFlag)
//			{
//				// Note: last parameter was not used
//				//                phiV = vortexPhi(PN,Al,a,b,s,local.row(0),local.row(1),local.row(2));
//				phiV = vortexPhi(PN, Al, a, b, s, local.row(0), local.row(1));
//			}
//			phi += srcSidePhi(PN, Al, phiV, a, b, s);
//		}
//		phi /= (4 * M_PI);
//	}
//}

void bodyPanel::panelVInf(const Eigen::Vector3d &POI, Eigen::Vector3d &vSrc,Eigen::Vector3d &vDub)
{
    
    // VSAero source and doublet velocity influence formulation
    Eigen::Vector3d pjk = POI-center;
    Eigen::Matrix3d local = getLocalSys();
    double PN = pjk.dot(local.row(2));
    if (pjk.norm()/longSide > 5)
    {
        vSrc = -pntSrcV(pjk); // Connor added this negative. Figure 10.17 (Katz) shows a positive velocity influence from pnt source. CPanel matches this in tests, however, sign convention for CPanel/VSaero is opposite of Katz
        vDub = pntDubV(local.row(2),pjk);
    }
    else
    {
        Eigen::Vector3d p1,p2,a,b,s,l,m,n;
        l = local.row(0);
        m = local.row(1);
        n = local.row(2);
        pjk = POI-center;
        double Al;
        for (nodes_index_type i=0; i<nodes.size(); i++)
        {
            if (i!=nodes.size()-1)
            {
                p1 = nodes[i]->getPnt();
                p2 = nodes[i+1]->getPnt();
            }
            else
            {
                p1 = nodes[i]->getPnt();
                p2 = nodes[0]->getPnt();
            }
            a = POI-p1;
            b = POI-p2;
            s = a-b;
            Al = n.dot(s.cross(a));
            
            vDub += vortexV(a,b,s);
            vSrc += srcSideV(PN,Al,a,b,s,l,m,n);
        }
        vDub /= (4*M_PI);
        vSrc /= (4*M_PI);
    }
}

Eigen::Vector3d bodyPanel::pntVInf(const Eigen::Vector3d &POI){
    // Function should only be used for CPanel test function.
    
    // VSAero source and doublet velocity influence formulation
    Eigen::Vector3d pjk = POI-center;
    Eigen::Matrix3d local = getLocalSys();
    
    return sourceStrength*pntSrcV(pjk) + doubletStrength*pntDubV(local.row(2),pjk); //pntDubVinf
}


double bodyPanel::srcSidePhi(const double &PN,const double &Al, const double &phiV,const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s)
{
    double A,B,S;
    A = a.norm();
    B = b.norm();
    S = s.norm();
    double GL = 0;
    if (std::abs(A+B-S) > 0 && S > 0)
    {
    	GL = 1/S*log(std::abs((A+B+S)/(A+B-S)));
    }
    return (Al*GL-PN*phiV);
}

Eigen::Vector3d bodyPanel::srcSideV(const double &PN,const double &Al, const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s,const Eigen::Vector3d &l,const Eigen::Vector3d &m,const Eigen::Vector3d &n)
{
    double A,B,S;
    A = a.norm();
    B = b.norm();
    S = s.norm();
    double GL = 0;
    if (std::abs(A+B-S) > 0 && S > 0)
    {
        GL = 1/S*log(std::abs((A+B+S)/(A+B-S)));
    }
// NOTE: last parameter is not used
//    double CJK = vortexPhi(PN,Al,a,b,s,l,m,n);
    double CJK = vortexPhi(PN,Al,a,b,s,l,m);
    return (GL*(s.dot(m)*l-s.dot(l)*m)+CJK*n);
}

inline double bodyPanel::pntSrcPhi(const double &PJK)
{
    return area/(4*M_PI*PJK);
}

inline Eigen::Vector3d bodyPanel::pntSrcV(const Eigen::Vector3d &pjk)
{
    return area*pjk/(4*M_PI*pow(pjk.norm(),3));
}

void bodyPanel::setCluster()
{
    int dim;
    size_t buffer = 10;
    if (tipFlag)
    {
        dim = 2;
    }
    else
    {
        dim = 3;
    }
    
    size_t nObs = static_cast<size_t>(chtlsnd::factorial(TSorder+dim)/(chtlsnd::factorial(dim)*chtlsnd::factorial(TSorder))-1)+buffer; // Binomial Coefficient
    size_t nPanels = (nObs+1)/2; // ceil(nObs/2) for integers
    bodyPanels_index_type oldSize = cluster.size();
    cluster.push_back(this);
    bool upFlag = upper;
    bool lowFlag = lower;
    while (cluster.size() < nPanels)
    {
        std::vector<bodyPanel*> toAdd;
        for (bodyPanels_index_type i=oldSize; i<cluster.size(); i++)
        {
            std::vector<bodyPanel*> temp = cluster[i]->getNeighbors();
            for (bodyPanels_index_type j=0; j<temp.size(); j++)
            {
                if (clusterTest(temp[j], 5*M_PI/6,upFlag,lowFlag))
                {
                    // Do not include panels on other side of discontinuity (i.e. wake), panels already in cluster, or this panel
                    if (temp[j]->isUpper())
                    {
                        upFlag = true;
                    }
                    else if (temp[j]->isLower())
                    {
                        lowFlag = true;
                    }
                    if (std::find(toAdd.begin(), toAdd.end(), temp[j]) == toAdd.end())
                    {
                        toAdd.push_back(temp[j]);
                    }
                }
            }
        }
        
        if (toAdd.size() == 0 && cluster.size() < nPanels)
        {
            // No valid panels could be found so decrease order of taylor series
            while (cluster.size() < nPanels+1 && TSorder >= 1)
            {
                if (buffer > 0)
                {
                    buffer -= 1;
                }
                else
                {
                    TSorder -= 1;
                }
                
                // Do Something to fix error when cluster can't be found.
                
                nObs = static_cast<size_t>(chtlsnd::factorial(TSorder+dim)/(chtlsnd::factorial(dim)*chtlsnd::factorial(TSorder))) + buffer; // Binomial Coefficient
                nPanels = (nObs+1)/2; // ceil(nObs/2) for integers
                
            }
            
            if (cluster.size() > 0) {
                cluster.erase(cluster.begin()+static_cast<bodyPanels_type::difference_type>(nPanels+1),cluster.end());
            }
            break;
        }
        else
        {
            oldSize = cluster.size();
            for (bodyPanels_index_type i=0; i<toAdd.size(); i++)
            {
                cluster.push_back(toAdd[i]);
                if (cluster.size() == nPanels+1)
                {
                    break;
                }
            }
        }
    }
}


Eigen::Vector3d bodyPanel::pntVelocity(const Eigen::Vector3d &pnt,double pntPotential, double PG, const Eigen::Vector3d &Vinf)
{
    Eigen::Vector3d vel;
    
    if (cluster.size() == 0)
    {
        setCluster();
    }
    
    if (tipFlag)
    {
        vel = velocity2D(pnt,pntPotential);
    }
    else
    {
        vel = velocity3D(pnt,pntPotential);
        if (vel != vel)
        {
            // NaN can be returned if supporting data was all on a flat patch (essentially 2D)
            vel = velocity2D(pnt,pntPotential);
        }
    }
    
    if (vel != vel)
    {
        // Last resort approximation if CHTLS not able to compute velocity. From Kinney (CBAERO).
        vel = (bezNormal.cross(Vinf)).cross(bezNormal);
    }
    
    vel(0) /= PG; // Prandlt Glauert Correction
    return vel;
}


Eigen::Vector3d bodyPanel::velocity2D(const Eigen::Vector3d &pnt,double pntPotential)
{
    int dim = 2;
    std::vector<bodyPanel*> clust = cluster;
    if (pnt == center)
    {
        clust.erase(clust.begin());
    }
    Eigen::Vector3d vel;
    Eigen::MatrixXd Xf,Xb,Vb;
    Eigen::VectorXd df;
    
    Eigen::Vector3d pntLocal = global2local(pnt,true);
    Eigen::Vector2d X0 = pntLocal.head(2);
    Eigen::Vector2d V0 = Eigen::Vector2d::Zero();
    Xf.resize(static_cast<Eigen::MatrixXd::Index>(clust.size()),3);
    Xb = Eigen::MatrixXd::Zero(0,dim);
    Vb = Eigen::MatrixXd::Zero(0,dim);
    df.resize(static_cast<Eigen::MatrixXd::Index>(clust.size()));
    for (bodyPanels_index_type i=0; i<clust.size(); i++)
    {
        Xf.row(static_cast<Eigen::MatrixXd::Index>(i)) = global2local(clust[i]->getCenter(),true);
        df(static_cast<Eigen::MatrixXd::Index>(i)) = clust[i]->getPotential()-pntPotential;
    }
    Eigen::MatrixXd xLocal = Xf.block(0,0,static_cast<Eigen::MatrixXd::Index>(clust.size()),2);
    chtlsnd tipV(X0,xLocal,TSorder,Xb,Vb,V0);
    Eigen::Vector3d vLocal;
    vLocal(0) = tipV.getF().row(0)*df;
    vLocal(1) = tipV.getF().row(1)*df;
    vLocal(2) = 0;
    vel = local2global(vLocal,false);
    
    return vel;
}

Eigen::Vector3d bodyPanel::velocity3D(const Eigen::Vector3d &pnt,double pntPotential)
{
    int dim = 3;
    std::vector<bodyPanel*> clust = cluster;
    if (pnt == center)
    {
        clust.erase(clust.begin());
    }
    Eigen::Vector3d vel;
    Eigen::MatrixXd Xf,Xb,Vb;
    Eigen::VectorXd df;
    
    Xf.resize(static_cast<Eigen::MatrixXd::Index>(clust.size()),dim);
    Xb.resize(static_cast<Eigen::MatrixXd::Index>(clust.size()),dim);
    Vb.resize(static_cast<Eigen::MatrixXd::Index>(clust.size()),dim);
    df.resize(static_cast<Eigen::MatrixXd::Index>(clust.size()));
    for (bodyPanels_index_type i=0; i<clust.size(); i++)
    {
        Xf.row(static_cast<Eigen::MatrixXd::Index>(i)) = clust[i]->getCenter();
        Vb.row(static_cast<Eigen::MatrixXd::Index>(i)) = clust[i]->getBezNormal();
        df(static_cast<Eigen::MatrixXd::Index>(i)) = clust[i]->getPotential()-pntPotential;
    }
    Xb = Xf;
    chtlsnd vWeights(pnt,Xf,TSorder,Xb,Vb,Eigen::Vector3d::Zero());
    
    vel(0) = vWeights.getF().row(0)*df;
    vel(1) = vWeights.getF().row(1)*df;
    vel(2) = vWeights.getF().row(2)*df;
    
    return vel;
}

void bodyPanel::computeCp(double Vinf)
{
    Cp = (1-pow(velocity.norm()/Vinf,2));
}

void bodyPanel::computeCp(double Vinf,double dt){
    double dPhi_dt = ( prevPotential - potential ) / dt;
    
    Cp = 1 - pow( velocity.norm()/Vinf , 2) - 2/(Vinf*Vinf) * dPhi_dt; // Katz 13.168
}

void bodyPanel::computeVelocity(double PG, const Eigen::Vector3d &Vinf)
{
    velocity = pntVelocity(center,potential,PG,Vinf);
	//std::cout << velocity << "\n" << std::endl;
}



Eigen::Vector3d bodyPanel::computeMoments(const Eigen::Vector3d &cg)
{
    Eigen::Vector3d r = center-cg;
    Eigen::Vector3d F = -Cp*bezNormal*area;
    return r.cross(F);
}

bool bodyPanel::clusterTest(bodyPanel* other,double angle,bool upFlag,bool lowFlag)
{
    if ((upFlag && other->isLower()) || (lowFlag && other->isUpper()))
    {
        return false;
    }
    
    if (tipFlag != other->isTipPan())
    {
        return false;
    }
    
    double dot = other->getNormal().dot(normal);
    // Floating point error can cause panels with the same normal vector to result in a dot product greater than one, causing acos to return nan.
    if (dot > 1)
    {
        dot = 1;
    }
    else if (dot < -1)
    {
        dot = -1;
    }
    return (acos(dot) < angle && std::find(cluster.begin(),cluster.end(),other)==cluster.end() && other != this);
}

void bodyPanel::setLSflag()
{
    parentSurf->setLSflag();
}

bool bodyPanel::nearTrailingEdge()
{
    // Returns true if panel is on trailing edge or a neighbor is on trailing edge
    if (parentSurf->isLiftingSurf())
    {
        if (upper || lower)
        {
            return true;
        }
        else
        {
            std::vector<bodyPanel*> neighbs = getNeighbors();
            for (bodyPanels_index_type i=0; i<neighbs.size(); i++)
            {
                if (neighbs[i]->isUpper() || neighbs[i]->isLower())
                {
                    return true;
                }
            }
        }
    }
    return false;
}

double bodyPanel::dist2Pan(bodyPanel* other)
{
    return (other->getCenter()-center).norm();
}

std::vector<bodyPanel*> bodyPanel::getRelatedPanels()
{
    std::vector<bodyPanel*> pans;
    pans.push_back(this);
    std::vector<bodyPanel*> nodePans;
    
    for (nodes_index_type i=0; i<nodes.size(); i++)
    {
        nodePans = nodes[i]->getBodyPans();
        for (size_t j=0; j<nodePans.size(); j++)
        {
            if (std::find(pans.begin(),pans.end(),nodePans[j]) == pans.end())
            {
                pans.push_back(nodePans[j]);
            }
            
        }
    }
    return pans;
}

Eigen::Vector3d bodyPanel::partStretching(particle* part){
    
    Eigen::Vector3d partStretching;
    double dist2panel = (part->pos-center).norm();
    Eigen::MatrixXd velGradMat = Eigen::Matrix3d::Zero();
    
    // Far field convention usually either 3 or 5, see Katz or Chris' thesis
    bool isFarField = false;
    if(dist2panel/longSide > 5){
        isFarField = true;
    }
    
    if(isFarField)
    {
        velGradMat += velocityGradientPointSource(part->pos);
        velGradMat += velocityGradientPointDoublet(part->pos);
    }
    else if (this->getNodes().size() == 3)
    {
        velGradMat += velocityGradientTriSource(part->pos);
        velGradMat += velocityGradientDoublet(part->pos);
    }
    else
    {
        velGradMat += velocityGradientQuadSource(part->pos);
        velGradMat += velocityGradientDoublet(part->pos);
    }
    
    partStretching = velGradMat*part->strength;
    
    return partStretching;
}



Eigen::Matrix3d bodyPanel::velocityGradientPointSource(Eigen::Vector3d POI){
    // Equations in documentation
    Eigen::Matrix3d velGradMat;
    // dudx  dvdx  dwdx
    // dudy  dvdy  dwdy
    // dudz  dvdz  dwdz
    
    double x, y, z, x0, y0, z0;
    
    x = this->getCenter().x(); y = this->getCenter().y(); z = this->getCenter().z();
    
    x0 = POI.x(); y0 = POI.y(); z0 = POI.z();
    
    double sdd = pow(pow(x-x0,2) + pow(y-y0,2) + z*z,2.5); // source deriv denom
    
    double derConst = sourceStrength*area/(4*M_PI);
    velGradMat(0,0) = (-2*(x-x0)*(x-x0) + (y-y0)*(y-y0) + z*z)/sdd;
    velGradMat(1,0) = -3*(x-x0)*(y-y0)/sdd;
    velGradMat(2,0) = -3*(x-x0)*z/sdd;
    
    velGradMat(0,1) = -3*(x-x0)*(y-y0)/sdd;
    velGradMat(1,1) = ((x-x0)*(x-x0)- (y-y0)*(y-y0) + z*z)/sdd;
    velGradMat(2,1) = -3*z*(y-y0)/sdd;
    
    velGradMat(0,2) = -3*(x-x0)*(z-z0)/sdd;
    velGradMat(1,2) = -3*(y-y0)*(z-z0)/sdd;
    velGradMat(2,2) = (-3*z*(z-z0) + ((x-x0)*(x-x0) + (y-y0)*(y-y0) + z*z))/sdd;
    
    velGradMat *= derConst;
    
    return velGradMat;
}


Eigen::Matrix3d bodyPanel::velocityGradientQuadSource(Eigen::Vector3d POI){
    Eigen::Matrix3d velGradMat;
    // dudx  dvdx  dwdx
    // dudy  dvdy  dwdy
    // dudz  dvdz  dwdz
    
    Eigen::Vector3d n1global = this->getNodes()[0]->getPnt();
    Eigen::Vector3d n2global = this->getNodes()[1]->getPnt();
    Eigen::Vector3d n3global = this->getNodes()[2]->getPnt();
    Eigen::Vector3d n4global = this->getNodes()[3]->getPnt();
    
    Eigen::Vector3d n1 = global2local(n1global, true);
    Eigen::Vector3d n2 = global2local(n2global, true);
    Eigen::Vector3d n3 = global2local(n3global, true);
    Eigen::Vector3d n4 = global2local(n4global, true);
    Eigen::Vector3d POIloc = global2local(POI, true);
    
    double x1, y1, x2, y2, x3, y3, x4, y4;
    x1 = n1.x(); x2 = n2.x(); x3 = n3.x(), x4 = n4.x();
    y1 = n1.y(); y2 = n2.y(); y3 = n3.y(), y4 = n4.y();
    
    double x, y, z;
    x = POIloc.x(); y = POIloc.y(); z = POIloc.z();
    
    
    double derConst = sourceStrength/(4*M_PI);
    
    // Terms
    double d12,d23,d34,d41;
    d12 = pow((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1),0.5);
    d23 = pow((x3-x2)*(x3-x2) + (y3-y2)*(y3-y2),0.5);
    d34 = pow((x4-x3)*(x4-x3) + (y4-y3)*(y4-y3),0.5);
    d41 = pow((x1-x4)*(x1-x4) + (y1-y4)*(y1-y4),0.5);
    
    double r1,r2,r3,r4;
    r1 = pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
    r2 = pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
    r3 = pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
    r4 = pow((x-x4)*(x-x4) + (y-y4)*(y-y4) + z*z,0.5);
    
    double dr1dx, dr2dx, dr3dx, dr4dx;
    dr1dx = (x-x1)/pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
    dr2dx = (x-x2)/pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
    dr3dx = (x-x3)/pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
    dr4dx = (x-x4)/pow((x-x4)*(x-x4) + (y-y4)*(y-y4) + z*z,0.5);
    
    double dr1dy, dr2dy, dr3dy, dr4dy;
    dr1dy = (y-y1)/pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
    dr2dy = (y-y2)/pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
    dr3dy = (y-y3)/pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
    dr4dy = (y-y4)/pow((x-x4)*(x-x4) + (y-y4)*(y-y4) + z*z,0.5);
    
    double dr1dz, dr2dz, dr3dz, dr4dz;
    dr1dz = (z)/pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
    dr2dz = (z)/pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
    dr3dz = (z)/pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
    dr4dz = (z)/pow((x-x4)*(x-x4) + (y-y4)*(y-y4) + z*z,0.5);
    
    double m12,m23,m34,m41;
    m12 = (y2-y1)/(x2-x1);
    m23 = (y3-y2)/(x3-x2);
    m34 = (y4-y3)/(x4-x3);
    m41 = (y1-y4)/(x1-x4);
    
    double e1,e2,e3,e4;
    e1 = (x-x1)*(x-x1) + z*z;
    e2 = (x-x2)*(x-x2) + z*z;
    e3 = (x-x3)*(x-x3) + z*z;
    e4 = (x-x4)*(x-x4) + z*z;
    
    double h1,h2,h3,h4;
    h1 = (x-x1)*(y-y1);
    h2 = (x-x2)*(y-y2);
    h3 = (x-x3)*(y-y3);
    h4 = (x-x4)*(y-y4);
    
    double de1dx,de2dx,de3dx,de4dx;
    de1dx = 2*(x-x1);
    de2dx = 2*(x-x2);
    de3dx = 2*(x-x3);
    de4dx = 2*(x-x4);
    
    double dh1dx,dh2dx,dh3dx,dh4dx;
    dh1dx = (y-y1);
    dh2dx = (y-y2);
    dh3dx = (y-y3);
    dh4dx = (y-y4);
    
    double dh1dy,dh2dy,dh3dy,dh4dy;
    dh1dy = (x-x1);
    dh2dy = (x-x2);
    dh3dy = (x-x3);
    dh4dy = (x-x4);
    
    double de1dz,de2dz,de3dz,de4dz;
    de1dz = 2*z;
    de2dz = 2*z;
    de3dz = 2*z;
    de4dz = 2*z;
    
    
    velGradMat(0,0)=(y2-y1)*(r1+r2+d12)*(dr1dx + dr2dx)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (y3-y2)*(r2+r3+d23)*(dr2dx + dr3dx)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (y4-y3)*(r3+r4+d34)*(dr3dx + dr4dx)*(2*d34)/(d34)/(r3+r4-d34)/(pow(r3+r4+d34,2))+
    (y1-y4)*(r4+r1+d41)*(dr4dx + dr1dx)*(2*d41)/(d41)/(r4+r1-d41)/(pow(r4+r1+d41,2));
    
    velGradMat(1,0)=(y2-y1)*(r1+r2+d12)*(dr1dy + dr2dy)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (y3-y2)*(r2+r3+d23)*(dr2dy + dr3dy)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (y4-y3)*(r3+r4+d34)*(dr3dy + dr4dy)*(2*d34)/(d34)/(r3+r4-d34)/(pow(r3+r4+d34,2))+
    (y1-y4)*(r4+r1+d41)*(dr4dy + dr1dy)*(2*d41)/(d41)/(r4+r1-d41)/(pow(r4+r1+d41,2));
    
    velGradMat(2,0)=(y2-y1)*(r1+r2+d12)*(dr1dz + dr2dz)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (y3-y2)*(r2+r3+d23)*(dr2dz + dr3dz)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (y4-y3)*(r3+r4+d34)*(dr3dz + dr4dz)*(2*d34)/(d34)/(r3+r4-d34)/(pow(r3+r4+d34,2))+
    (y1-y4)*(r4+r1+d41)*(dr4dz + dr1dz)*(2*d41)/(d41)/(r4+r1-d41)/(pow(r4+r1+d41,2));
    
    velGradMat(0,1)=(x1-x2)*(r1+r2+d12)*(dr1dx + dr2dx)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (x2-x3)*(r2+r3+d23)*(dr2dx + dr3dx)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (x3-x4)*(r3+r4+d34)*(dr3dx + dr4dx)*(2*d34)/(d34)/(r3+r4-d34)/(pow(r3+r4+d34,2))+
    (x4-x1)*(r4+r1+d41)*(dr4dx + dr1dx)*(2*d41)/(d41)/(r4+r1-d41)/(pow(r4+r1+d41,2));
    
    velGradMat(1,1)=(x1-x2)*(r1+r2+d12)*(dr1dy + dr2dy)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (x2-x3)*(r2+r3+d23)*(dr2dy + dr3dy)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (x3-x4)*(r3+r4+d34)*(dr3dy + dr4dy)*(2*d34)/(d34)/(r3+r4-d34)/(pow(r3+r4+d34,2))+
    (x4-x1)*(r4+r1+d41)*(dr4dy + dr1dy)*(2*d41)/(d41)/(r4+r1-d41)/(pow(r4+r1+d41,2));
    
    velGradMat(2,1)=(x1-x2)*(r1+r2+d12)*(dr1dz + dr2dz)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (x2-x3)*(r2+r3+d23)*(dr2dz + dr3dz)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (x3-x4)*(r3+r4+d34)*(dr3dz + dr4dz)*(2*d34)/(d34)/(r3+r4-d34)/(pow(r3+r4+d34,2))+
    (x4-x1)*(r4+r1+d41)*(dr4dz + dr1dz)*(2*d41)/(d41)/(r4+r1-d41)/(pow(r4+r1+d41,2));
    
    velGradMat(0,2)=(1/(1+pow((m12*e1-h1)/(z*r1),2)))*((z*r1*(m12*de1dx - dh1dx) - (m12*e1-h1)*(z*dr1dx))/pow(z*r1,2))-
    (1/(1+pow((m12*e2-h2)/(z*r2),2)))*((z*r2*(m12*de2dx - dh2dx) - (m12*e2-h2)*(z*dr2dx))/pow(z*r2,2))+
    (1/(1+pow((m23*e2-h2)/(z*r2),2)))*((z*r2*(m23*de2dx - dh2dx) - (m23*e2-h2)*(z*dr2dx))/pow(z*r2,2))-
    (1/(1+pow((m23*e3-h3)/(z*r3),2)))*((z*r3*(m23*de3dx - dh3dx) - (m23*e3-h3)*(z*dr3dx))/pow(z*r3,2))+
    (1/(1+pow((m34*e3-h3)/(z*r3),2)))*((z*r3*(m34*de3dx - dh3dx) - (m34*e3-h3)*(z*dr3dx))/pow(z*r3,2))-
    (1/(1+pow((m34*e4-h4)/(z*r4),2)))*((z*r4*(m34*de4dx - dh4dx) - (m34*e4-h4)*(z*dr4dx))/pow(z*r4,2))+
    (1/(1+pow((m41*e4-h4)/(z*r4),2)))*((z*r4*(m41*de4dx - dh4dx) - (m41*e4-h4)*(z*dr4dx))/pow(z*r4,2))-
    (1/(1+pow((m41*e1-h1)/(z*r1),2)))*((z*r1*(m41*de1dx - dh1dx) - (m41*e1-h1)*(z*dr1dx))/pow(z*r1,2));
    
    velGradMat(1,2)=(1/(1+pow((m12*e1-h1)/(z*r1),2)))*((-z*r1*dh1dy - (m12*e1-h1)*(z*dr1dy))/pow(z*r1,2))-
    (1/(1+pow((m12*e2-h2)/(z*r2),2)))*((-z*r2*dh2dy - (m12*e2-h2)*(z*dr2dy))/pow(z*r2,2))+
    (1/(1+pow((m23*e2-h2)/(z*r2),2)))*((-z*r2*dh2dy - (m23*e2-h2)*(z*dr2dy))/pow(z*r2,2))-
    (1/(1+pow((m23*e3-h3)/(z*r3),2)))*((-z*r3*dh3dy - (m23*e3-h3)*(z*dr3dy))/pow(z*r3,2))+
    (1/(1+pow((m34*e3-h3)/(z*r3),2)))*((-z*r3*dh3dy - (m34*e3-h3)*(z*dr3dy))/pow(z*r3,2))-
    (1/(1+pow((m34*e4-h4)/(z*r4),2)))*((-z*r4*dh4dy - (m34*e4-h4)*(z*dr4dy))/pow(z*r4,2))+
    (1/(1+pow((m41*e4-h4)/(z*r4),2)))*((-z*r4*dh4dy - (m41*e4-h4)*(z*dr4dy))/pow(z*r4,2))-
    (1/(1+pow((m41*e1-h1)/(z*r1),2)))*((-z*r1*dh1dy - (m41*e1-h1)*(z*dr1dy))/pow(z*r1,2));
    
    velGradMat(2,2)=(1/(1+pow((m12*e1-h1)/(z*r1),2)))*((z*r1*m12*de1dz - (m12*e1-h1)*(r1+z*dr1dz))/pow(z*r1,2))-
    (1/(1+pow((m12*e2-h2)/(z*r2),2)))*((z*r2*m12*de2dz - (m12*e2-h2)*(r2+z*dr2dz))/pow(z*r2,2))+
    (1/(1+pow((m23*e2-h2)/(z*r2),2)))*((z*r2*m23*de2dz - (m23*e2-h2)*(r2+z*dr2dz))/pow(z*r2,2))-
    (1/(1+pow((m23*e3-h3)/(z*r3),2)))*((z*r3*m23*de3dz - (m23*e3-h3)*(r3+z*dr3dz))/pow(z*r3,2))+
    (1/(1+pow((m34*e3-h3)/(z*r3),2)))*((z*r3*m34*de3dz - (m34*e3-h3)*(r3+z*dr3dz))/pow(z*r3,2))-
    (1/(1+pow((m34*e4-h4)/(z*r4),2)))*((z*r4*m34*de4dz - (m34*e4-h4)*(r4+z*dr4dz))/pow(z*r4,2))+
    (1/(1+pow((m41*e4-h4)/(z*r4),2)))*((z*r4*m41*de4dz - (m41*e4-h4)*(r4+z*dr4dz))/pow(z*r4,2))-
    (1/(1+pow((m41*e1-h1)/(z*r1),2)))*((z*r1*m41*de1dz - (m41*e1-h1)*(r1+z*dr1dz))/pow(z*r1,2));
    
    velGradMat *= derConst;
    
    return velGradMat;
}

Eigen::Matrix3d bodyPanel::velocityGradientTriSource(Eigen::Vector3d POI){
    Eigen::Matrix3d velGradMat;
    // dudx  dvdx  dwdx
    // dudy  dvdy  dwdy
    // dudz  dvdz  dwdz
    
    Eigen::Vector3d n1global = this->getNodes()[0]->getPnt();
    Eigen::Vector3d n2global = this->getNodes()[1]->getPnt();
    Eigen::Vector3d n3global = this->getNodes()[2]->getPnt();
    
    Eigen::Vector3d n1 = global2local(n1global, true);
    Eigen::Vector3d n2 = global2local(n2global, true);
    Eigen::Vector3d n3 = global2local(n3global, true);
    Eigen::Vector3d POIloc = global2local(POI, true);
    
    double x1, y1, x2, y2, x3, y3;
    x1 = n1.x(); x2 = n2.x(); x3 = n3.x();
    y1 = n1.y(); y2 = n2.y(); y3 = n3.y();
    
    double x, y, z;
    x = POIloc.x(); y = POIloc.y(); z = POIloc.z();
    
    double derConst = sourceStrength/(4*M_PI);
    
    // Terms
    double d12,d23,d31;
    d12 = pow((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1),0.5);
    d23 = pow((x3-x2)*(x3-x2) + (y3-y2)*(y3-y2),0.5);
    d31 = pow((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3),0.5);
    
    double r1,r2,r3;
    r1 = pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
    r2 = pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
    r3 = pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
    
    double dr1dx, dr2dx, dr3dx;
    dr1dx = (x-x1)/pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
    dr2dx = (x-x2)/pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
    dr3dx = (x-x3)/pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
    
    double dr1dy, dr2dy, dr3dy;
    dr1dy = (y-y1)/pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
    dr2dy = (y-y2)/pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
    dr3dy = (y-y3)/pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
    
    double dr1dz, dr2dz, dr3dz;
    dr1dz = (z)/pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z,0.5);
    dr2dz = (z)/pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z,0.5);
    dr3dz = (z)/pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + z*z,0.5);
    
    double m12,m23,m31;
    m12 = (y2-y1)/(x2-x1);
    m23 = (y3-y2)/(x3-x2);
    m31 = (y1-y3)/(x1-x3);
    
    double e1,e2,e3;
    e1 = (x-x1)*(x-x1) + z*z;
    e2 = (x-x2)*(x-x2) + z*z;
    e3 = (x-x3)*(x-x3) + z*z;
    
    double h1,h2,h3;
    h1 = (x-x1)*(y-y1);
    h2 = (x-x2)*(y-y2);
    h3 = (x-x3)*(y-y3);
    
    double de1dx,de2dx,de3dx;
    de1dx = 2*(x-x1);
    de2dx = 2*(x-x2);
    de3dx = 2*(x-x3);
    
    double dh1dx,dh2dx,dh3dx;
    dh1dx = (y-y1);
    dh2dx = (y-y2);
    dh3dx = (y-y3);
    
    double dh1dy,dh2dy,dh3dy;
    dh1dy = (x-x1);
    dh2dy = (x-x2);
    dh3dy = (x-x3);
    
    double de1dz,de2dz,de3dz;
    de1dz = 2*z;
    de2dz = 2*z;
    de3dz = 2*z;
    
    
    velGradMat(0,0)=(y2-y1)*(r1+r2+d12)*(dr1dx + dr2dx)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (y3-y2)*(r2+r3+d23)*(dr2dx + dr3dx)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (y1-y3)*(r3+r1+d31)*(dr3dx + dr1dx)*(2*d31)/(d31)/(r3+r1-d31)/(pow(r3+r1+d31,2));
    
    velGradMat(1,0)=(y2-y1)*(r1+r2+d12)*(dr1dy + dr2dy)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (y3-y2)*(r2+r3+d23)*(dr2dy + dr3dy)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (y1-y3)*(r3+r1+d31)*(dr3dy + dr1dy)*(2*d31)/(d31)/(r3+r1-d31)/(pow(r3+r1+d31,2));
    
    velGradMat(2,0)=(y2-y1)*(r1+r2+d12)*(dr1dz + dr2dz)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (y3-y2)*(r2+r3+d23)*(dr2dz + dr3dz)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (y1-y3)*(r3+r1+d31)*(dr3dz + dr1dz)*(2*d31)/(d31)/(r3+r1-d31)/(pow(r3+r1+d31,2));
    
    velGradMat(0,1)=(x1-x2)*(r1+r2+d12)*(dr1dx + dr2dx)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (x2-x3)*(r2+r3+d23)*(dr2dx + dr3dx)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (x3-x1)*(r3+r1+d31)*(dr3dx + dr1dx)*(2*d31)/(d31)/(r3+r1-d31)/(pow(r3+r1+d31,2));
    
    velGradMat(1,1)=(x1-x2)*(r1+r2+d12)*(dr1dy + dr2dy)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (x2-x3)*(r2+r3+d23)*(dr2dy + dr3dy)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (x3-x1)*(r3+r1+d31)*(dr3dy + dr1dy)*(2*d31)/(d31)/(r3+r1-d31)/(pow(r3+r1+d31,2));
    
    velGradMat(2,1)=(x1-x2)*(r1+r2+d12)*(dr1dz + dr2dz)*(2*d12)/(d12)/(r1+r2-d12)/(pow(r1+r2+d12,2))+
    (x2-x3)*(r2+r3+d23)*(dr2dz + dr3dz)*(2*d23)/(d23)/(r2+r3-d23)/(pow(r2+r3+d23,2))+
    (x3-x1)*(r3+r1+d31)*(dr3dz + dr1dz)*(2*d31)/(d31)/(r3+r1-d31)/(pow(r3+r1+d31,2));
    
    velGradMat(0,2)=(1/(1+pow((m12*e1-h1)/(z*r1),2)))*((z*r1*(m12*de1dx - dh1dx) - (m12*e1-h1)*(z*dr1dx))/pow(z*r1,2))-
    (1/(1+pow((m12*e2-h2)/(z*r2),2)))*((z*r2*(m12*de2dx - dh2dx) - (m12*e2-h2)*(z*dr2dx))/pow(z*r2,2))+
    (1/(1+pow((m23*e2-h2)/(z*r2),2)))*((z*r2*(m23*de2dx - dh2dx) - (m23*e2-h2)*(z*dr2dx))/pow(z*r2,2))-
    (1/(1+pow((m23*e3-h3)/(z*r3),2)))*((z*r3*(m23*de3dx - dh3dx) - (m23*e3-h3)*(z*dr3dx))/pow(z*r3,2))+
    (1/(1+pow((m31*e3-h3)/(z*r3),2)))*((z*r3*(m31*de3dx - dh3dx) - (m31*e3-h3)*(z*dr3dx))/pow(z*r3,2))-
    (1/(1+pow((m31*e1-h1)/(z*r1),2)))*((z*r1*(m31*de1dx - dh1dx) - (m31*e1-h1)*(z*dr1dx))/pow(z*r1,2));
    
    velGradMat(1,2)=(1/(1+pow((m12*e1-h1)/(z*r1),2)))*((-z*r1*dh1dy - (m12*e1-h1)*(z*dr1dy))/pow(z*r1,2))-
    (1/(1+pow((m12*e2-h2)/(z*r2),2)))*((-z*r2*dh2dy - (m12*e2-h2)*(z*dr2dy))/pow(z*r2,2))+
    (1/(1+pow((m23*e2-h2)/(z*r2),2)))*((-z*r2*dh2dy - (m23*e2-h2)*(z*dr2dy))/pow(z*r2,2))-
    (1/(1+pow((m23*e3-h3)/(z*r3),2)))*((-z*r3*dh3dy - (m23*e3-h3)*(z*dr3dy))/pow(z*r3,2))+
    (1/(1+pow((m31*e3-h3)/(z*r3),2)))*((-z*r3*dh3dy - (m31*e3-h3)*(z*dr3dy))/pow(z*r3,2))-
    (1/(1+pow((m31*e1-h1)/(z*r1),2)))*((-z*r1*dh1dy - (m31*e1-h1)*(z*dr1dy))/pow(z*r1,2));
    
    velGradMat(2,2)=(1/(1+pow((m12*e1-h1)/(z*r1),2)))*((z*r1*m12*de1dz - (m12*e1-h1)*(r1+z*dr1dz))/pow(z*r1,2))-
    (1/(1+pow((m12*e2-h2)/(z*r2),2)))*((z*r2*m12*de2dz - (m12*e2-h2)*(r2+z*dr2dz))/pow(z*r2,2))+
    (1/(1+pow((m23*e2-h2)/(z*r2),2)))*((z*r2*m23*de2dz - (m23*e2-h2)*(r2+z*dr2dz))/pow(z*r2,2))-
    (1/(1+pow((m23*e3-h3)/(z*r3),2)))*((z*r3*m23*de3dz - (m23*e3-h3)*(r3+z*dr3dz))/pow(z*r3,2))+
    (1/(1+pow((m31*e3-h3)/(z*r3),2)))*((z*r3*m31*de3dz - (m31*e3-h3)*(r3+z*dr3dz))/pow(z*r3,2))-
    (1/(1+pow((m31*e1-h1)/(z*r1),2)))*((z*r1*m31*de1dz - (m31*e1-h1)*(r1+z*dr1dz))/pow(z*r1,2));
    
    velGradMat *= derConst;
    
    return velGradMat;
}


void bodyPanel::linComputeVelocity(double PG, Eigen::Vector3d &Vinf)
{
	double mu_x, mu_y;
	Eigen::Vector3d vertDubStrengths, linDubConsts, panVel;
	Eigen::Matrix3d vertsMat = linVertsMatrix(true);
	vertDubStrengths = linGetDubStrengths();

	linDubConsts = vertsMat.inverse() * vertDubStrengths;
	mu_x = linDubConsts[1];
	mu_y = linDubConsts[2];

	panVel[0] = mu_x;
	panVel[1] = mu_y;
	panVel[2] = 0;

	velocity = local2global(panVel, false);
	velocity(0) /= PG;
}


Eigen::Vector3d bodyPanel::linComputeVelocity2(const double PG, Eigen::Vector3d Vinf)
{
	double s, beta2, denom, a;
	Eigen::Vector3d vertsPhiAvg, coeffsPhiAvg, gradPhiAvg, gradDubs;
	Eigen::Vector3d vAvgTan, vAvgNorm, vDiffTan, vDiffNorm, vAvg, vDiff;
	Eigen::Vector3d pertVel, velLower, myPanPertVel, myPertVel;
	Eigen::Matrix3d vertsMat;

	vertsMat = linVertsMatrix(true);
	//vertsMat = supVertsMatrix(supLocalNodes);

	// vertex based average potential
	vertsPhiAvg = linGetDubStrengths();

	// linear panel equation average potential
	coeffsPhiAvg = vertsMat.inverse() * vertsPhiAvg;

	// gradient of average potential
	gradPhiAvg << coeffsPhiAvg[1], coeffsPhiAvg[2], 0;

	// gradient of 'raw' doublet strength coeffs
	gradDubs << linDubCoeffs[1], linDubCoeffs[2], 0;

	// Tangent and Normal components of vAvg and vDiff
	vAvgTan = supTransMat.transpose() * gradPhiAvg;
	vAvgNorm = (0.5*sourceStrength / denom) * normal;

	vDiffTan = supTransMat.transpose() * gradDubs;
	vDiffNorm = (sourceStrength / denom) * normal;

	// Final average velocity and velocity difference
	vAvg = vAvgTan + vAvgNorm;
	vDiff = vDiffTan + vDiffNorm;

	velLower = vAvg - vDiff / 2; // should be 0!

								 // Compute perturbation velocity
	pertVel = vAvg + vDiff / 2; // upper velocity (i.e. perturbation velocity)

								// Compute total velocity
	velocity = Vinf + pertVel;

	// Return perturbation velocity for use in Cp calcs
	return pertVel;
}


Eigen::Vector3d bodyPanel::supComputeVelocity(Eigen::Vector3d Vinf, const double mach, bool velCorrection)
{
	double s, beta2, denom, a;
	Eigen::Vector3d vertsPhiAvg, coeffsPhiAvg, gradPhiAvg, gradDubs;
	Eigen::Vector3d vAvgTan, vAvgNorm, vDiffTan, vDiffNorm, vAvg, vDiff;
	Eigen::Vector3d pertVel, velLower, myPanPertVel, myPertVel;
	Eigen::Matrix3d vertsMat, B;

	vertsMat = supVertsMatrix(supLocalNodes);

	// vertex based average potential
	vertsPhiAvg = linGetDubStrengths();

	// linear panel equation average potential
	coeffsPhiAvg = vertsMat.inverse() * vertsPhiAvg;

	// gradient of average potential
	gradPhiAvg << coeffsPhiAvg[1], coeffsPhiAvg[2], 0;

	// gradient of 'raw' doublet strength coeffs
	gradDubs << linDubCoeffs[1], linDubCoeffs[2], 0;

	s = -1; // -1 for supersonic, +1 for subsonic
	beta2 = s * (1 - pow(mach, 2));
	B << s*beta2, 0, 0, 0, 1, 0, 0, 0, 1;
	denom = normal.dot(B * normal); // normal dotted with conormal

	// Tangent and Normal components of vAvg and vDiff
	vAvgTan = supTransMat.transpose() * gradPhiAvg;
	vAvgNorm = (0.5*sourceStrength / denom) * normal;

	vDiffTan = supTransMat.transpose() * gradDubs;
	vDiffNorm = (sourceStrength / denom) * normal;

	// Final average velocity and velocity difference
	vAvg = vAvgTan + vAvgNorm;
	vDiff = vDiffTan + vDiffNorm;

	velLower = vAvg - vDiff / 2; // should be 0!

	// Compute perturbation velocity
	pertVel = vAvg + vDiff / 2; // upper velocity (i.e. perturbation velocity)
	
	// Compute total velocity
	velocity = Vinf + pertVel;

	a = Vinf.norm() / mach;
	panelMach = velocity.norm() / a;

	// Return perturbation velocity for use in Cp calcs
	return pertVel;
}


void bodyPanel::supComputeCp(Eigen::Vector3d Vinf, const double mach, Eigen::Vector3d pertVelWind, Eigen::Vector3d pertVelBody)
{
	// 1st Order Cp (Linearized)
	Cp1 = -2 * pertVelWind.x() / Vinf.norm();

	// 2nd Order Cp (Slender Body Assumption)
	Cp2s = (-2 * pertVelWind.x() / Vinf.norm()) - ((pow(pertVelWind.y(), 2) + pow(pertVelWind.z(), 2)) / pow(Vinf.norm(), 2));

	// 2nd Order Cp (no Slender Body Assumption)
	Cp2 = (-2 * pertVelWind.x() / Vinf.norm()) - (((1-pow(mach,2))*pow(pertVelWind.x(),2) + pow(pertVelWind.y(), 2) + pow(pertVelWind.z(), 2)) / pow(Vinf.norm(), 2));

	// Slender approximation in body axes for force calcs
	Cp = (-2 * pertVelBody.x() / Vinf.norm()) - ((pow(pertVelBody.y(), 2) + pow(pertVelBody.z(), 2)) / pow(Vinf.norm(), 2));
}

void bodyPanel::supFixCp(const double CpIn)
{
	Cp = CpIn;
}

void bodyPanel::supFixCp2s(const double CpIn)
{
	Cp2s = CpIn;
}

void bodyPanel::supFixMach(const double machIn)
{
	panelMach = machIn;
}



void bodyPanel::supFlipNormal()
{
	normal *= -1.0;
}


bool bodyPanel::supDODcheck(Eigen::Vector3d &P, const double Mach, Eigen::Vector3d &windDir)
{
	// P is POI, Q is panel center

	bool DODflag = false;
	Eigen::Vector3d Q, Pc, Qc, triDiams;
	double B, RBcntr, x0, y0, dist, triRad;
	B = sqrt(pow(Mach,2) - 1);
	Q = center;

	// Hyperbolic distance from POI to panel center
	RBcntr = pow(P.x() - Q.x(), 2) - pow(B, 2)*pow(P.y() - Q.y(), 2) - pow(B, 2)*pow(P.z() - Q.z(), 2);

	// Check the distance from the Mach cone boundary to the panel. If near, look more closely. If far away, check whether inside or out

	// CSYS shift for panel-to-Mach cone distance check
	Pc = P - P;
	Qc = Q - P;
	x0 = (Qc - Pc).dot(windDir);
	y0 = sqrt(pow((Qc - Pc).norm(), 2) - pow(x0, 2));
	if (y0 >= B * x0)
	{
		dist = sqrt(pow(x0 + B * y0, 2) / (1 + pow(B, 2)));
	}
	else
	{
		dist = sqrt(pow((Qc - Pc).norm(), 2));
	}

	for (size_t i = 0; i < nodes.size(); i++)
	{
		if (i != nodes.size() - 1)
		{
			triDiams[i] = (nodes[i]->getPnt() - nodes[i + 1]->getPnt()).norm();
		}
		else
		{
			triDiams[i] = (nodes[0]->getPnt() - nodes[i]->getPnt()).norm();
		}
	}
	triRad = triDiams.maxCoeff() / 2; // Panel radius, checked against panel distance from Mach cone

	if ((Pc - Qc).dot(windDir) >= 0 || dist < triRad) // Check if POI is downstream of panel, or if it's close to the panel
	{
		if (dist > triRad) // check if P is further from Mach cone than panel radius
		{
			if (RBcntr >= 0)
			{
				return DODflag = true;
			}
		}
		else // Panel is close to POI, look closer
		{
			double RBvert;
			Eigen::Vector3d vert;
			for (size_t i = 0; i < nodes.size(); i++) // Check if any panel vertices are inside the Mach cone
			{
				vert = nodes[i]->getPnt();
				RBvert = pow(P.x() - vert.x(), 2) - pow(B, 2)*pow(P.y() - vert.y(), 2) - pow(B, 2)*pow(P.z() - vert.z(), 2);
				if (RBvert >= 0 && (P - vert).dot(windDir) >= 0)
				{
					return DODflag = true;
				}
			}

			// If this point is reached, panel is either just outside of the Mach cone, or it intersects the Mach cone with no vertices inside the Mach cone
			// Thus, a closer check is needed. This is done in supPhiInf
			return DODflag = true;
		}
	}

	return DODflag;
}


void bodyPanel::supTransformPanel(double alpha, double beta, const double mach)
{
	// alpha and beta are input in degrees the transformed to radians

	alpha *= M_PI / 180;
	beta *= M_PI / 180;
	supSetG2LSmatrix(alpha, beta, mach);

	for (size_t i = 0; i < nodes.size(); i++)
	{
		supLocalNodes.push_back(supTransMat * (nodes[i]->getPnt() - center));
	}
}


void bodyPanel::linTransformPanel()
{
	for (size_t i = 0; i < nodes.size(); i++)
	{
		linLocalNodes.push_back(global2local(nodes[i]->getPnt(), true));
	}
}


 void bodyPanel::supSetG2LSmatrix(const double a, const double b, const double mach)
{
	 // Transformation from global CSYS to local-scaled CSYS

	// alpha and beta are input in radians
	Eigen::Vector3d cRef, nRef, nWind, scaleV, nWindScaled, vRef, uRef, col0, col1, col2;
	Eigen::Matrix3d ref2wind, machScaleB, machScaleC, machScaleWind;
	double beta, s, r;

	if (a != 0 || b != 0)
	{
		// Freestream direction in reference CSYS
		cRef << cos(a)*cos(b), -sin(b), sin(a)*cos(b);

		// Reference to wind transformation
		ref2wind << cos(a)*cos(b), -sin(b), sin(a)*cos(b),
			cos(a)*sin(b), cos(b), sin(a)*sin(b),
			-sin(a), 0, cos(a);
	}
	else
	{
		cRef << 1, 0, 0;
		ref2wind.setIdentity();
	}

	// Matrices to remove B (i.e. Mach) from integration
	s = -1.0;
	beta = sqrt(s * (1 - pow(mach, 2)));
	machScaleB << s * pow(beta, 2), 0, 0,
		0, 1, 0,
		0, 0, 1;
	machScaleC << 1, 0, 0,
		0, s*pow(beta, 2), 0,
		0, 0, s*pow(beta, 2);

	// Combined wind and Mach scaling matrix
	machScaleWind = ref2wind.transpose() * machScaleC * ref2wind;

	// Panel normal in reference CSYS
	nRef = normal;

	// Stuff
	vRef = nRef.cross(cRef) / (nRef.cross(cRef)).norm();
	uRef = vRef.cross(nRef);
	r = 1.0;
	if (signbit(nRef.dot(machScaleB * nRef)))
	{
		r = -1.0;
	}

	col0 = (1.0 / sqrt(abs(nRef.dot(machScaleB * nRef)))) * machScaleWind * uRef;
	col1 = (r * s / beta) * machScaleWind * vRef;
	col2 = (beta * nRef) / sqrt(abs(nRef.dot(machScaleB * nRef)));

	supTransMat.col(0) = col0;
	supTransMat.col(1) = col1;
	supTransMat.col(2) = col2;

	supTransMat.transposeInPlace();

	// Panel area correction term applied to computed source coefficients
	supAreaCorrect = 1.0 / (beta * sqrt(abs(nRef.dot(machScaleB * nRef))));
}


void bodyPanel::supPhiInf(const Eigen::Vector3d &POI, Eigen::Matrix<double, 1, Eigen::Dynamic> &Arow, double &srcPhi, bool DODflag, const double Mach)
{
	//if (DODflag)
	if (DODflag && !supInclinedFlag) // maybe make different supinclined panel func to make this more efficient
	{
		// Panel transformation projects original panel onto x-y plane w/ origin at panel center
		// Thus panel z-coords are zero, and freestream direction is [1 0 0]

		// Initialize
		double x, y, z, x1, y1, x2, y2;
		double m, lam, xm, xmc, ym1, ym2, ym1c, ym2c, s1, s2, sm1, sm2, r1, r2, R1, R2;
		double Bmach, Q1sum, srcCoeff;
		Eigen::Vector2d integralCoeffs; // [Q1, w0]
		Eigen::Vector3d localPOI, localWindDir, dubCoeffs, dubCoeffMat, dubVertsPhi;
		Eigen::Matrix3d dubVertsMat;
		bool mFlag, zFlag;

		localWindDir << 1, 0, 0;
		zFlag = false;
		Bmach = sqrt(pow(Mach, 2) - 1);
		integralCoeffs.setZero();
		dubCoeffs.setZero();
		srcCoeff = 0;

		double epsGenStrong, epsGenWeak;
		epsGenStrong = 1.0e-10;
		epsGenWeak = 1.0e-6;

		// Transform control point
		localPOI = supTransMat * (POI - center);
		x = localPOI.x();
		y = localPOI.y();
		z = localPOI.z();

		// Iterate through panel edges
		for (nodes_index_type i = 0; i < nodes.size(); i++)
		{
			integralCoeffs.setZero();
			mFlag = false;
			Eigen::Vector3d p1;
			Eigen::Vector3d p2;
			double eps1; // set based on offset of control point from it's node
			double eps2; // same as above
			if (i != nodes.size() - 1)
			{
				p1 = supLocalNodes[i];
				p2 = supLocalNodes[i + 1];

				eps1 = 2.0*nodes[i]->linGetCPoffset();
				eps2 = 2.0*nodes[i + 1]->linGetCPoffset();
			}
			else
			{
				p1 = supLocalNodes[i];
				p2 = supLocalNodes[0];

				eps1 = 2.0*nodes[i]->linGetCPoffset();
				eps2 = 2.0*nodes[0]->linGetCPoffset();
			}

			x1 = p1.x();
			y1 = p1.y();
			x2 = p2.x();
			y2 = p2.y();

			// Compute inf. coeff. geometric quantities
			m = (y2 - y1) / (x2 - x1);
			s1 = y - y1;
			s2 = y - y2;
			r1 = sqrt(pow(s1, 2) + pow(z, 2));
			r2 = sqrt(pow(s2, 2) + pow(z, 2));

			if (abs(m) < epsGenStrong) // Edge parallel to freestream
			{
				mFlag = true;
				m = 0;
				xmc = -s1;
				ym1c = -(x - x1);
				ym2c = -(x - x2);
				sm1 = x - x1;
				sm2 = x - x2;

				xm = xmc;
			}
			else if (abs(m) > 1.0 / epsGenStrong) // Edge perinducular to freestream
			{
				lam = 0;
				xm = (x - x1) - (y - y1)*lam;
				sm1 = xm + s1 * lam;
				sm2 = xm + s2 * lam;
				ym1 = s1 - sm1 * lam;
				ym2 = s2 - sm2 * lam;
			}
			else
			{
				lam = 1 / m;
				xm = (x - x1) - (y - y1) / m;
				sm1 = xm + s1 / m;
				sm2 = xm + s2 / m;
				ym1 = s1 - sm1 / m;
				ym2 = s2 - sm2 / m;
				xmc = m * xm;
				ym1c = m * ym1;
				ym2c = m * ym2;
			}

			// Check whether edge endpoints are inside or outside Mach cone. If outside (R has imaginary part), set to 0
			R1 = 0;
			R2 = 0;
			if ((localPOI - p1).dot(localWindDir) >= 0)
			{
				R1 = sqrt(pow(sm1, 2) - pow(r1, 2));
				if (isnan(R1))
				{
					R1 = 0;
				}
			}
			if ((localPOI - p2).dot(localWindDir) >= 0)
			{
				R2 = sqrt(pow(sm2, 2) - pow(r2, 2));
				if (isnan(R2))
				{
					R2 = 0;
				}
			}

			// Check that the edge is either inside or intersects the Mach cone
			if (R1 > abs(eps1) || R2 > abs(eps2))
			{
				if (abs(z) > epsGenStrong)
				{
					if (abs(1 - abs(m)) < epsGenWeak) // sonic edge (parallel to Mach cone)
					{
						//supOutputGeom(POI, true);
						integralCoeffs = supEdgeInfSon(ym1, ym2, xm, R1, R2, lam, z, eps1, eps2);
					}
					else if (abs(m) < 1) // subsonic edge (edge slope is less than that of the Mach cone)
					//if (abs(m) < 1) // subsonic edge (edge slope is less than that of the Mach cone)
					{
						integralCoeffs = supEdgeInfSub(R1, R2, ym1c, ym2c, xmc, m, z, eps1, eps2, mFlag);
					}
					else if (abs(m) > 1) // supersonic edge (edge slope is greater than that of the Mach cone)
					{
						integralCoeffs = supEdgeInfSup(R1, R2, ym1, ym2, xm, lam, z, eps1, eps2);
					}
				}
				else // z = 0
				{
					//--------------------------------------------------------------------------------------//
					// Need to be tested more thoroughly
					zFlag = true;
					if (abs(1 - abs(m)) < epsGenWeak) // sonic edge
					{
						double stuff = 0;
					}
					else if (abs(m) < 1) // subsonic edge
					{
						integralCoeffs = supEdgeInfSubInPlane(R1, R2, ym1c, ym2c, xmc, m, z, eps1, eps2, mFlag);
					}
					else if (abs(m) > 1) // supersonic edge
					{
						integralCoeffs = supEdgeInfSupInPlane(R1, R2, ym1, ym2, xm, lam, z, eps1, eps2);
					}
					//--------------------------------------------------------------------------------------//
				}
			} // Edge either intersects Mach cone with no edge points inside (i.e. inside Mach wedge), or is completely outside. Former case can only occur with supersonic edge
			else if (abs(m) > 1)
			{
				// Check that edge is upstream from POI, POI is between the endpoints of the edge, and it's within the Mach cone in the z-direction
				double wedgeCheck = pow(xm, 2) + (pow(z, 2) / pow(m, 2)) - pow(z, 2);
				if ((xm > 0) && (signbit(ym1) != signbit(ym2)) && (wedgeCheck >= 0))
				{
					double Q11, w01, Q12, w02, Q1, w0;

					double s1, t1;
					s1 = 1.0;
					t1 = 1.0;
					if (signbit(z * ym1))
					{
						s1 = -1.0;
					}
					if (signbit(ym1))
					{
						t1 = -1.0;
					}

					Q11 = s1 * M_PI / 2.0;
					w01 = t1 * M_PI / (2.0 * sqrt(1.0 - pow(lam, 2)));

					double s2, t2;
					s2 = 1.0;
					t2 = 1.0;
					if (signbit(z * ym2))
					{
						s2 = -1.0;
					}
					if (signbit(ym2))
					{
						t2 = -1.0;
					}
					Q12 = s2 * M_PI / 2.0;
					w02 = t2 * M_PI / (2.0 * sqrt(1.0 - pow(lam, 2)));

					Q1 = Q12 - Q11;
					w0 = w02 - w01;
					integralCoeffs[0] = Q1;
					integralCoeffs[1] = w0;
				}
			}

			// Sum doublet coefficients contributions from each edge
				// dubCoeffs = [Q1 w0 w0/m]
			dubCoeffs[0] += integralCoeffs[0];
			dubCoeffs[1] += integralCoeffs[1];
			if (!mFlag)
			{
				dubCoeffs[2] += integralCoeffs[1] / m;
			}

			// Sum source coefficient contributions from each edge
				// srcCoeff = [xm*w0]
			srcCoeff += xm * integralCoeffs[1];
		}

		double srcNum;
		if (zFlag)
		{
			dubCoeffs[1] = 0;
			dubCoeffs[2] = 0;
			srcNum = srcCoeff;
		}
		else
		{
			dubCoeffs[1] = -z * dubCoeffs[1];
			dubCoeffs[2] = -z * dubCoeffs[2];
			Q1sum = dubCoeffs[0];
			srcNum = srcCoeff - z * Q1sum;
		}

		srcPhi = (srcNum / (2.0 * M_PI)) * supAreaCorrect;

		dubCoeffs[0] = -dubCoeffs[0];
		dubCoeffs = dubCoeffs / (2.0 * M_PI);

		// Build matrix of doublet coefficients
		dubCoeffMat[0] = dubCoeffs[0];
		dubCoeffMat[1] = x * dubCoeffs[0] - dubCoeffs[1];
		dubCoeffMat[2] = y * dubCoeffs[0] - dubCoeffs[2];

		// Convert linear doublet equation to vertex based linear doublet equations
		dubVertsMat = supVertsMatrix(supLocalNodes);

		// Compute influence of each vertex on the field point
		dubVertsPhi = dubCoeffMat.transpose() * dubVertsMat.inverse();

		// Fill A matrix
		// Need to make this more efficient
		for (size_t i = 0; i < nodes.size(); i++)
		{
			// Check if influencing node is the same as the influenced node
			if (abs((POI - nodes[i]->getPnt()).norm()) < 2.0*nodes[i]->linGetCPoffset())
			{
				Arow[nodes[i]->getIndex()] = -0.5;
			}
			else
			{
				//--------------------------------------------------------------------------------------//
				// This is for subsonic LE, flat bottom cases
				/*if (normal.z() >= 0)
				{
					Arow[nodes[i]->getIndex()] += dubVertsPhi(i);
				}*/
				//--------------------------------------------------------------------------------------//
				Arow[nodes[i]->getIndex()] += dubVertsPhi(i);
			}
		}
	}
	else
	{
		for (size_t i = 0; i < nodes.size(); i++)
		{
			// Check if influencing node is the same as the influenced node
			if (abs((POI - nodes[i]->getPnt()).norm()) < 2.0*nodes[i]->linGetCPoffset())
			{
				Arow[nodes[i]->getIndex()] = -0.5;
			}
		}
	}
}



void bodyPanel::supPhiInf11012019(const Eigen::Vector3d &POI, Eigen::Matrix<double, 1, Eigen::Dynamic> &Arow, double &srcPhi, bool DODflag, const double Mach, size_t ii)
{
	//if (DODflag)
	if (DODflag && !supInclinedFlag) // maybe make different supinclined panel func to make this more efficient
	{
		// Panel transformation projects original panel onto x-y plane w/ origin at panel center
		// Thus panel z-coords are zero, and freestream direction is [1 0 0]

		// Initialize
		double x, y, z, x1, y1, x2, y2;
		double m, lam, xm, xmc, ym1, ym2, ym1c, ym2c, s1, s2, sm1, sm2, r1, r2, R1, R2;
		double Bmach, Q1sum, srcCoeff;
		Eigen::Vector2d integralCoeffs; // [Q1, w0]
		Eigen::Vector3d localPOI, localWindDir, dubCoeffs, dubCoeffMat, dubVertsPhi;
		Eigen::Matrix3d dubVertsMat;
		bool mFlag, zFlag;

		localWindDir << 1, 0, 0;
		zFlag = false;
		Bmach = sqrt(pow(Mach, 2) - 1);
		integralCoeffs.setZero();
		dubCoeffs.setZero();
		srcCoeff = 0;

		double epsGenStrong, epsGenWeak;
		epsGenStrong = 1.0e-10;
		epsGenWeak = 1.0e-6;

		// Transform control point
		localPOI = supTransMat * (POI - center);
		x = localPOI.x();
		y = localPOI.y();
		z = localPOI.z();

		// Iterate through panel edges
		for (nodes_index_type i = 0; i < nodes.size(); i++)
		{
			integralCoeffs.setZero();
			mFlag = false;
			Eigen::Vector3d p1;
			Eigen::Vector3d p2;
			double eps1; // set based on offset of control point from it's node
			double eps2; // same as above
			if (i != nodes.size() - 1)
			{
				p1 = supLocalNodes[i];
				p2 = supLocalNodes[i + 1];

				eps1 = 2.0*nodes[i]->linGetCPoffset();
				eps2 = 2.0*nodes[i + 1]->linGetCPoffset();
			}
			else
			{
				p1 = supLocalNodes[i];
				p2 = supLocalNodes[0];

				eps1 = 2.0*nodes[i]->linGetCPoffset();
				eps2 = 2.0*nodes[0]->linGetCPoffset();
			}

			x1 = p1.x();
			y1 = p1.y();
			x2 = p2.x();
			y2 = p2.y();

			// Compute inf. coeff. geometric quantities
			m = (y2 - y1) / (x2 - x1);
			s1 = y - y1;
			s2 = y - y2;
			r1 = sqrt(pow(s1, 2) + pow(z, 2));
			r2 = sqrt(pow(s2, 2) + pow(z, 2));

			if (abs(m) < epsGenStrong) // Edge parallel to freestream
			{
				mFlag = true;
				m = 0;
				xmc = -s1;
				ym1c = -(x - x1);
				ym2c = -(x - x2);
				sm1 = x - x1;
				sm2 = x - x2;

				xm = xmc;
			}
			else if (abs(m) > 1.0 / epsGenStrong) // Edge perinducular to freestream
			{
				lam = 0;
				xm = (x - x1) - (y - y1)*lam;
				sm1 = xm + s1 * lam;
				sm2 = xm + s2 * lam;
				ym1 = s1 - sm1 * lam;
				ym2 = s2 - sm2 * lam;
			}
			else
			{
				lam = 1 / m;
				xm = (x - x1) - (y - y1) / m;
				sm1 = xm + s1 / m;
				sm2 = xm + s2 / m;
				ym1 = s1 - sm1 / m;
				ym2 = s2 - sm2 / m;
				xmc = m * xm;
				ym1c = m * ym1;
				ym2c = m * ym2;
			}

			// Check whether edge endpoints are inside or outside Mach cone. If outside (R has imaginary part), set to 0
			R1 = 0;
			R2 = 0;
			if ((localPOI - p1).dot(localWindDir) >= 0)
			{
				R1 = sqrt(pow(sm1, 2) - pow(r1, 2));
				if (isnan(R1))
				{
					R1 = 0;
				}
			}
			if ((localPOI - p2).dot(localWindDir) >= 0)
			{
				R2 = sqrt(pow(sm2, 2) - pow(r2, 2));
				if (isnan(R2))
				{
					R2 = 0;
				}
			}

			// Check that the edge is either inside or intersects the Mach cone
			if (R1 > abs(eps1) || R2 > abs(eps2))
			{
				if (abs(z) > epsGenStrong)
				{
					if (abs(1 - abs(m)) < epsGenWeak) // sonic edge (parallel to Mach cone)
					{
						//supOutputGeom(POI, true);
						integralCoeffs = supEdgeInfSon(ym1, ym2, xm, R1, R2, lam, z, eps1, eps2);
					}
					else if (abs(m) < 1) // subsonic edge (edge slope is less than that of the Mach cone)
										 //if (abs(m) < 1) // subsonic edge (edge slope is less than that of the Mach cone)
					{
						integralCoeffs = supEdgeInfSub(R1, R2, ym1c, ym2c, xmc, m, z, eps1, eps2, mFlag);
					}
					else if (abs(m) > 1) // supersonic edge (edge slope is greater than that of the Mach cone)
					{
						integralCoeffs = supEdgeInfSup(R1, R2, ym1, ym2, xm, lam, z, eps1, eps2);
					}
				}
				else // z = 0
				{
					//--------------------------------------------------------------------------------------//
					// Need to be tested more thoroughly
					zFlag = true;
					if (abs(1 - abs(m)) < epsGenWeak) // sonic edge
					{
						double stuff = 0;
					}
					else if (abs(m) < 1) // subsonic edge
					{
						integralCoeffs = supEdgeInfSubInPlane(R1, R2, ym1c, ym2c, xmc, m, z, eps1, eps2, mFlag);
					}
					else if (abs(m) > 1) // supersonic edge
					{
						integralCoeffs = supEdgeInfSupInPlane(R1, R2, ym1, ym2, xm, lam, z, eps1, eps2);
					}
					//--------------------------------------------------------------------------------------//
				}
			} // Edge either intersects Mach cone with no edge points inside (i.e. inside Mach wedge), or is completely outside. Former case can only occur with supersonic edge
			else if (abs(m) > 1)
			{
				// Check that edge is upstream from POI, POI is between the endpoints of the edge, and it's within the Mach cone in the z-direction
				double wedgeCheck = pow(xm, 2) + (pow(z, 2) / pow(m, 2)) - pow(z, 2);
				if ((xm > 0) && (signbit(ym1) != signbit(ym2)) && (wedgeCheck >= 0))
				{
					double Q11, w01, Q12, w02, Q1, w0;

					double s1, t1;
					s1 = 1.0;
					t1 = 1.0;
					if (signbit(z * ym1))
					{
						s1 = -1.0;
					}
					if (signbit(ym1))
					{
						t1 = -1.0;
					}

					Q11 = s1 * M_PI / 2.0;
					w01 = t1 * M_PI / (2.0 * sqrt(1.0 - pow(lam, 2)));

					double s2, t2;
					s2 = 1.0;
					t2 = 1.0;
					if (signbit(z * ym2))
					{
						s2 = -1.0;
					}
					if (signbit(ym2))
					{
						t2 = -1.0;
					}
					Q12 = s2 * M_PI / 2.0;
					w02 = t2 * M_PI / (2.0 * sqrt(1.0 - pow(lam, 2)));

					Q1 = Q12 - Q11;
					w0 = w02 - w01;
					integralCoeffs[0] = Q1;
					integralCoeffs[1] = w0;
				}
			}

			// Sum doublet coefficients contributions from each edge
			// dubCoeffs = [Q1 w0 w0/m]
			dubCoeffs[0] += integralCoeffs[0];
			dubCoeffs[1] += integralCoeffs[1];
			if (!mFlag)
			{
				dubCoeffs[2] += integralCoeffs[1] / m;
			}

			// Sum source coefficient contributions from each edge
			// srcCoeff = [xm*w0]
			srcCoeff += xm * integralCoeffs[1];
		}

		double srcNum;
		if (zFlag)
		{
			dubCoeffs[1] = 0;
			dubCoeffs[2] = 0;
			srcNum = srcCoeff;
		}
		else
		{
			dubCoeffs[1] = -z * dubCoeffs[1];
			dubCoeffs[2] = -z * dubCoeffs[2];
			Q1sum = dubCoeffs[0];
			srcNum = srcCoeff - z * Q1sum;
		}

		srcPhi = (srcNum / (2.0 * M_PI)) * supAreaCorrect;

		dubCoeffs[0] = -dubCoeffs[0];
		dubCoeffs = dubCoeffs / (2.0 * M_PI);

		// Build matrix of doublet coefficients
		dubCoeffMat[0] = dubCoeffs[0];
		dubCoeffMat[1] = x * dubCoeffs[0] - dubCoeffs[1];
		dubCoeffMat[2] = y * dubCoeffs[0] - dubCoeffs[2];

		// Convert linear doublet equation to vertex based linear doublet equations
		dubVertsMat = supVertsMatrix(supLocalNodes);

		// Compute influence of each vertex on the field point
		dubVertsPhi = dubCoeffMat.transpose() * dubVertsMat.inverse();

		// Fill A matrix
		// Need to make this more efficient
		for (size_t i = 0; i < nodes.size(); i++)
		{
			// Check if influencing node is the same as the influenced node
			if (ii == nodes[i]->getIndex())
			//if (abs((POI - nodes[i]->getPnt()).norm()) < 2.0*nodes[i]->linGetCPoffset())
			{
				Arow[nodes[i]->getIndex()] = -0.5;
			}
			else
			{
				//--------------------------------------------------------------------------------------//
				// This is for subsonic LE, flat bottom cases
				/*if (normal.z() >= 0)
				{
				Arow[nodes[i]->getIndex()] += dubVertsPhi(i);
				}*/
				//--------------------------------------------------------------------------------------//
				//Arow[nodes[i]->getIndex()] += dubVertsPhi(i);
			}
		}
	}
	else
	{
		for (size_t i = 0; i < nodes.size(); i++)
		{
			// Check if influencing node is the same as the influenced node
			if (ii == nodes[i]->getIndex())
			//if (abs((POI - nodes[i]->getPnt()).norm()) < 2.0*nodes[i]->linGetCPoffset())
			{
				Arow[nodes[i]->getIndex()] = -0.5;
			}
		}
	}
}



Eigen::Vector2d bodyPanel::supEdgeInfSon(const double ym1, const double ym2, const double xm, const double R1, const double R2, const double lam, const double z, const double eps1, const double eps2)
{
	Eigen::Vector2d integralCoeffs;
	double Q1, w0, ym1Neg, ym2Neg, zr;
	ym1Neg = -ym1;
	ym2Neg = -ym2;

	zr = (ym1Neg*R2 - ym2Neg * R1) / (ym1Neg*ym2Neg + (1 - pow(lam, 2))*R1*R2);
	w0 = zr * (1 - ((1 - pow(lam, 2)*pow(zr, 2)) / 3) + ((pow(1 - pow(lam, 2), 2)*pow(zr, 4)) / 5) - ((pow(1 - pow(lam, 2), 3)*pow(zr, 6)) / 7));

	if (R1 > eps1 && R2 > eps2) // Both edge points are inside the Mach cone
	{
		Q1 = atan2(z * xm * (ym1Neg*R2 - ym2Neg * R1), pow(z, 2)*ym2Neg*ym1Neg + pow(xm, 2)*R1*R2);
	}
	else // One edge point is inside Mach cone
	{
		double Q11, Q12;
		if (R1 > eps1)
		{
			Q11 = atan2(z*ym1, xm*R1);
		}
		else
		{
			double s1 = 1.0;
			if (signbit(z * ym1))
			{
				s1 = -1.0;
			}

			Q11 = s1 * M_PI / 2;
		}
		if (R2 > eps2)
		{
			Q12 = atan2(z*ym2, xm*R2);
		}
		else
		{
			double s2 = 1.0;
			if (signbit(z * ym2))
			{
				s2 = -1.0;
			}
			Q12 = s2 * M_PI / 2;
		}
		Q1 = Q12 - Q11;
	}

	integralCoeffs[0] = Q1;
	integralCoeffs[1] = w0;

	return integralCoeffs;
}


Eigen::Vector2d bodyPanel::supEdgeInfSub(const double R1, const double R2, const double ym1c, const double ym2c, const double xmc, const double m, const double z, const double eps1, const double eps2, bool mFlag)
{
	Eigen::Vector2d integralCoeffs;
	double Q1, w0, ym1cNeg, ym2cNeg;
	ym1cNeg = -ym1c;
	ym2cNeg = -ym2c;

	if (R1 > eps1 && R2 > eps2) // Both edge points are inside the Mach cone
	{
		Q1 = atan2(z*xmc * (ym1cNeg*R2 - ym2cNeg * R1), pow(z, 2) * ym1cNeg * ym2cNeg + pow(xmc, 2)*R1*R2);
		if (mFlag) // edge parallel to freestream
		{
			w0 = log((ym2cNeg + R2) / (ym1cNeg + R1));
		}
		else
		{
			w0 = (m / sqrt(1 - pow(m, 2))) * log((ym2cNeg + R2 * sqrt(1 - pow(m, 2))) / (ym1cNeg + R1 * sqrt(1 - pow(m, 2))));
		}
	}
	else // One edge point is inside Mach cone
	{
		double Q11, Q12, w01, w02, s;
		s = 1.0;
		if (signbit(z))
		{
			s = -1.0;
		}

		if (mFlag) // edge parallel to freestream
		{
			if (R1 > eps1)
			{
				Q11 = s * atan2(xmc * R1, -abs(z) * ym1c);
				w01 = (1.0 / (2 * sqrt(1 - pow(m, 2)))) * log((-ym1c + R1 * sqrt(1 - pow(m, 2))) / (-ym1c - R1 * sqrt(1 - pow(m, 2))));
			}
			else
			{
				Q11 = 0;
				w01 = 0;
			}
			if (R2 > eps2)
			{
				Q12 = s * atan2(xmc * R2, -abs(z) * ym2c);
				w02 = (1.0 / (2 * sqrt(1 - pow(m, 2)))) * log((-ym2c + R2 * sqrt(1 - pow(m, 2))) / (-ym2c - R2 * sqrt(1 - pow(m, 2))));
			}
			else
			{
				Q12 = 0;
				w02 = 0;
			}
		}
		else
		{
			if (R1 > eps1)
			{
				Q11 = s * atan2(xmc * R1, -abs(z) * ym1c);
				w01 = (m / (2 * sqrt(1 - pow(m, 2)))) * log((-ym1c + R1 * sqrt(1 - pow(m, 2))) / (-ym1c - R1 * sqrt(1 - pow(m, 2))));
			}
			else
			{
				Q11 = 0;
				w01 = 0;
			}
			if (R2 > eps2)
			{
				Q12 = s * atan2(xmc * R2, -abs(z) * ym2c);
				w02 = (m / (2 * sqrt(1 - pow(m, 2)))) * log((-ym2c + R2 * sqrt(1 - pow(m, 2))) / (-ym2c - R2 * sqrt(1 - pow(m, 2))));
			}
			else
			{
				Q12 = 0;
				w02 = 0;
			}
		}

		Q1 = Q12 - Q11;
		w0 = w02 - w01;
	}

	integralCoeffs[0] = Q1;
	integralCoeffs[1] = w0;

	return integralCoeffs;
}


Eigen::Vector2d bodyPanel::supEdgeInfSup(const double R1, const double R2, const double ym1, const double ym2, const double xm, const double lam, const double z, const double eps1, const double eps2)
{
	Eigen::Vector2d integralCoeffs;
	double Q1, w0, ym1Neg, ym2Neg;
	ym1Neg = -ym1;
	ym2Neg = -ym2;

	if (R1 > eps1 && R2 > eps2) // Both edge points are inside the Mach cone
	{
		Q1 = atan2(z * xm * (ym1Neg*R2 - ym2Neg * R1), pow(z,2)*ym2Neg*ym1Neg + pow(xm,2)*R1*R2);
		w0 = (1 / sqrt(1 - pow(lam, 2))) * atan2(sqrt(1 - pow(lam, 2)) * (ym1Neg*R2 - ym2Neg * R1), ym1Neg*ym2Neg + (1 - pow(lam, 2))*R1*R2);
	}
	else // One edge point is inside Mach cone
	{
		double Q11, Q12, w01, w02;
		if (R1 > eps1)
		{
			Q11 = atan2(z*ym1, xm*R1);
			w01 = (1 / sqrt(1 - pow(lam, 2))) * atan2(ym1, R1*sqrt(1 - pow(lam, 2)));
		}
		else
		{
			double s1, t1;
			s1 = 1.0;
			t1 = 1.0;
			if (signbit(z * ym1))
			{
				s1 = -1.0;
			}
			if (signbit(ym1))
			{
				t1 = -1.0;
			}

			Q11 = s1 * M_PI / 2;
			w01 = t1 * M_PI / (2 * sqrt(1 - pow(lam, 2)));
		}
		if (R2 > eps2)
		{
			Q12 = atan2(z*ym2, xm*R2);
			w02 = (1 / sqrt(1 - pow(lam, 2))) * atan2(ym2, R2*sqrt(1 - pow(lam, 2)));
		}
		else
		{
			double s2, t2;
			s2 = 1.0;
			t2 = 1.0;
			if (signbit(z * ym2))
			{
				s2 = -1.0;
			}
			if (signbit(ym2))
			{
				t2 = -1.0;
			}
			Q12 = s2 * M_PI / 2;
			w02 = t2 * M_PI / (2 * sqrt(1 - pow(lam, 2)));
		}
		Q1 = Q12 - Q11;
		w0 = w02 - w01;
	}

	integralCoeffs[0] = Q1;
	integralCoeffs[1] = w0;

	return integralCoeffs;
}


Eigen::Vector2d bodyPanel::supEdgeInfSubInPlane(const double R1, const double R2, const double ym1c, const double ym2c, const double xmc, const double m, const double z, const double eps1, const double eps2, bool mFlag)
{
	Eigen::Vector2d integralCoeffs;
	double Q1, w0, ym1cNeg, ym2cNeg;
	ym1cNeg = -ym1c;
	ym2cNeg = -ym2c;

	if (R1 > eps1 && R2 > eps2) // Both edge points are inside the Mach cone
	{
		Q1 = 0;

		if (mFlag) // edge parallel to freestream
		{
			//w0 = log((ym2cNeg + R2) / (ym1cNeg + R1));
			w0 = 0;
		}
		else
		{
			w0 = (m / sqrt(1 - pow(m, 2))) * log((ym2cNeg + R2 * sqrt(1 - pow(m, 2))) / (ym1cNeg + R1 * sqrt(1 - pow(m, 2))));
		}
	}
	else // One edge point is inside Mach cone
	{
		double Q11, Q12, w01, w02, s;
		s = 1.0;
		if (signbit(z*xmc))
		{
			s = -1.0;
		}

		Q11 = 0;
		Q12 = 0;
		w01 = 0;
		w02 = 0;
		if (mFlag) // edge parallel to freestream
		{
			if (R1 > eps1)
			{
				Q11 = s * M_PI / 2;
				//w01 = (1.0 / (2 * sqrt(1 - pow(m, 2)))) * log((-ym1c + R1 * sqrt(1 - pow(m, 2))) / (-ym1c - R1 * sqrt(1 - pow(m, 2))));
				w01 = 0;
			}
			if (R2 > eps2)
			{
				Q12 = s * M_PI / 2;
				//w02 = (1.0 / (2 * sqrt(1 - pow(m, 2)))) * log((-ym2c + R2 * sqrt(1 - pow(m, 2))) / (-ym2c - R2 * sqrt(1 - pow(m, 2))));
				w02 = 0;
			}
		}
		else
		{
			if (R1 > eps1)
			{
				Q11 = s * M_PI / 2;
				//w01 = (m / (2 * sqrt(1 - pow(m, 2)))) * log((-ym1c + R1 * sqrt(1 - pow(m, 2))) / (-ym1c - R1 * sqrt(1 - pow(m, 2))));
				w01 = 0;
			}
			if (R2 > eps2)
			{
				Q12 = s * M_PI / 2;
				//w02 = (m / (2 * sqrt(1 - pow(m, 2)))) * log((-ym2c + R2 * sqrt(1 - pow(m, 2))) / (-ym2c - R2 * sqrt(1 - pow(m, 2))));
				w02 = 0;
			}
		}

		Q1 = Q12 - Q11;
		w0 = w02 - w01;

		if (isnan(Q1) || isnan(w0))
		{
			double stuff = 0;
		}
	}

	integralCoeffs[0] = Q1;
	integralCoeffs[1] = w0;

	return integralCoeffs;
}


Eigen::Vector2d bodyPanel::supEdgeInfSupInPlane(const double R1, const double R2, const double ym1, const double ym2, const double xm, const double lam, const double z, const double eps1, const double eps2)
{
	Eigen::Vector2d integralCoeffs;
	double Q1, w0, ym1Neg, ym2Neg;
	ym1Neg = -ym1;
	ym2Neg = -ym2;

	Q1 = 0;

	if (R1 > eps1 && R2 > eps2) // Both edge points are inside the Mach cone
	{
		w0 = (1 / sqrt(1 - pow(lam, 2))) * atan2(sqrt(1 - pow(lam, 2)) * (ym1Neg*R2 - ym2Neg * R1), ym1Neg*ym2Neg + (1 - pow(lam, 2))*R1*R2);
	}
	else // One edge point is inside Mach cone
	{
		double w01, w02;
		if (R1 > eps1)
		{
			w01 = (1 / sqrt(1 - pow(lam, 2))) * atan2(ym1, R1*sqrt(1 - pow(lam, 2)));
		}
		else
		{
			double t1 = 1.0;
			if (signbit(ym1))
			{
				t1 = -1.0;
			}

			w01 = t1 * M_PI / (2 * sqrt(1 - pow(lam, 2)));
		}
		if (R2 > eps2)
		{
			w02 = (1 / sqrt(1 - pow(lam, 2))) * atan2(ym2, R2*sqrt(1 - pow(lam, 2)));
		}
		else
		{
			double t2 = 1.0;
			if (signbit(ym2))
			{
				t2 = -1.0;
			}

			w02 = t2 * M_PI / (2 * sqrt(1 - pow(lam, 2)));
		}

		w0 = w02 - w01;
	}

	if (isnan(Q1) || isnan(w0))
	{
		double stuff = 0;
	}

	integralCoeffs[0] = Q1;
	integralCoeffs[1] = w0;

	return integralCoeffs;
}


void bodyPanel::supSetMu()
{
	Eigen::Vector3d vertDubStrengths;
	Eigen::Matrix3d vertsMat = supVertsMatrix(supLocalNodes);
	vertDubStrengths = linGetOrigDubStrengths();

	linDubCoeffs = vertsMat.inverse() * vertDubStrengths;
	// [mu_0 mu_x mu_y]
}


void bodyPanel::linSetMu()
{
	Eigen::Vector3d vertDubStrengths;
	Eigen::Matrix3d vertsMat = linVertsMatrix(true);
	vertDubStrengths = linGetDubStrengths();

	linDubCoeffs = vertsMat.inverse() * vertDubStrengths;
	// [mu_0 mu_x mu_y]
	std::cout << linDubCoeffs.y() << std::endl;
}


bool bodyPanel::supSuperinclinedCheck(const double B, Eigen::Matrix3d &body2wind)
{
	//bool isSupInclined = false;
	Eigen::Vector3d nWind, ncWind;

	nWind = body2wind * normal;
	ncWind = nWind;
	ncWind(0) *= -pow(B, 2);
	if (nWind.dot(ncWind) <= 0)
	{
		//isSupInclined = true;
		supInclinedFlag = true;
	}

	return supInclinedFlag;
}


///////////////////////////////////////////////////////////////////////////////// Funcion for testing against MATLAB
void bodyPanel::supOutputGeom(const Eigen::Vector3d &POI, bool outPOI)
{
	std::ofstream fid;
	std::string	myFile = "D:\\Desktop\\Thesis\\Code\\MATLAB\\Supersonic\\stuff.csv";
	fid.open(myFile);

	fid << nodes[0]->getPnt().x() << "," << nodes[0]->getPnt().y() << "," << nodes[0]->getPnt().z() << "\n";
	fid << nodes[1]->getPnt().x() << "," << nodes[1]->getPnt().y() << "," << nodes[1]->getPnt().z() << "\n";
	fid << nodes[2]->getPnt().x() << "," << nodes[2]->getPnt().y() << "," << nodes[2]->getPnt().z() << "\n";
	if (outPOI)
	{
		fid << POI.x() << "," << POI.y() << "," << POI.z();
		fid << "\n";
	}
	fid.close();
}


//Eigen::Vector3d bodyPanel::supVelCorrection(Eigen::Vector3d pertVel, const double mach)
//{
//	Eigen::Vector3d pertVelCorrected, pertMassFlux, totMassFlux;
//	double Bmach2, normLocDensity, lessCheck, moreCheck;
//	Bmach2 = pow(mach, 2) - 1; // Bmach2 = 1 w/ mach = sqrt(2)
//
//	pertMassFlux = pertVel;
//	pertMassFlux.x() *= Bmach2;
//	totMassFlux = velocity + pertMassFlux;
//
//	if (pertVel.x() < 0)
//	{
//		//velocity = totMassFlux / (1 - pow(mach, 2)*(pertVel.x()/Vinf.norm()));
//		pertVelCorrected = pertMassFlux / (1 - pow(mach, 2)*pertVel.x());
//	}
//
//	return pertVelCorrected;
//}
