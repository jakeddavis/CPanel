//
//  runCase.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 10/13/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "cpCase.h"

cpCase::~cpCase()
{
    for (int i=0; i<bStreamlines.size(); i++)
    {
        delete bStreamlines[i];
    }
}

void cpCase::run(bool printFlag, bool surfStreamFlag, bool stabDerivFlag, bool vortPartFlag)
{
    bool converged;
    std::string check = "\u2713";
    setSourceStrengths();
    
    converged = solveMatrixEq();

    std::cout << "starting wake" << std::endl;
    
    if(vortPartFlag){
        std::cout << "Writing timestep " << timeStep << " files..." << std::endl;
        
        
        compVelocity(); // This doesn't need to be here unless solving Cl/Cd etc
        writeFiles(); // Might take this out to or at least make it an option?
        timeStep++;
        
        
//        bool breakFlag = false;
//        std::cout << "time step = " << dt << std::endl;
        for(int i=0; i<numSteps; i++){
            std::cout << "Time step " << timeStep<< "/" << numSteps << ". Flow time = " << timeStep*dt << std::endl;

//            particleStrengthUpdate(); //commenting until big bug is found.
//            particleStrengthUpdateGaussian();
//            collapseBufferWake();
            collapseWakeForEachEdge();
            
            
            if(accelerate){
                
//                partOctree.clearTree();
//                partOctree.addData(particles);
                
                partOctree.removeData();
                partOctree.setMaxMembers(10); //Bares Hut
                partOctree.addData(particles);
                
                for (int i=0; i<particles.size(); i++) {
                    particles[i]->velOn = {0,0,0};
                }
                FMM.build(&partOctree);
            }
            
            std::cout << "Setting source strengths..." << std::endl;
            setSourceStrengths();
            std::cout << "Solving singularity strengths..." << std::endl;
            converged = solveVPmatrixEq();
//            std::cout << "Computing forces..." << std::endl;
//            compVelocity();
            
         

            std::cout << "writing files" << std::endl;
            writeFiles();
            timeStep++;
            

            
            
            std::cout << "Convecting " << particles.size() << " particles" << std::endl;
            convectParticles();

        }
        
//        std::cout << "New Trefftz Plane CD = " << trefftzPlaneCd(particles) << std::endl;
        
        std::cout << "CL=[";
        for(int i=0; i<CL.size(); i++){
            std::cout << CL[i] << ", ";
        }
        std::cout << "];" << std::endl;
    }
    
    if (printFlag)
    {
        std::cout << std::setw(17) << std::left << check << std::flush;
    }
    
    compVelocity();
    
    if (printFlag)
    {
        std::cout << std::setw(17) << std::left << check << std::flush;
    }
    
    trefftzPlaneAnalysis();
    
    if (printFlag)
    {
        std::cout << std::setw(18) << std::left << check << std::flush;
    }

    if (surfStreamFlag)
    {
        createStreamlines();
        if (printFlag)
        {
            std::cout << std::setw(16) << std::left << check << std::flush;
        }
    }
    else
    {
        if (printFlag)
        {
            std::cout << std::setw(16) << std::left << "X" << std::flush;
        }
    }
    
    if (stabDerivFlag)
    {
        stabilityDerivatives();
        if (printFlag)
        {
            std::cout << std::setw(23) << std::left << check << std::endl;
        }
    }
    else
    {
        if (printFlag)
        {
            std::cout << std::setw(23) << std::left << "X" << std::endl;
        }
    }
    
        
    if (!converged && printFlag)
    {
        std::cout << "*** Warning : Solution did not converge ***" << std::endl;;
    }
    
    if (printFlag){
        if(!vortPartFlag){ // Printed above
            writeFiles();
        }
    }
    
}

Eigen::Vector3d cpCase::windToBody(double V, double alpha, double beta)
{
    alpha *= M_PI/180;
    beta *= M_PI/180;
    
    Eigen::Matrix3d T;
    T(0,0) = cos(alpha)*cos(beta);
    T(0,1) = cos(alpha)*sin(beta);
    T(0,2) = -sin(alpha);
    T(1,0) = -sin(beta);
    T(1,1) = cos(beta);
    T(1,2) = 0;
    T(2,0) = sin(alpha)*cos(beta);
    T(2,1) = sin(alpha)*sin(beta);
    T(2,2) = cos(alpha);
    Eigen::Vector3d Vt;
    Eigen::Vector3d Vel;
    Vt << V,0,0;
    transform = T;

    Vel = transform*Vt;
    
    return Vel;
}

Eigen::Vector3d cpCase::bodyToWind(const Eigen::Vector3d &vec)
{
    Eigen::Matrix3d T = transform.transpose();
    return T*vec;
}

void cpCase::setSourceStrengths()
{
    sigmas.resize(bPanels->size());
    for (int i=0; i<bPanels->size(); i++)
    {
        Eigen::Vector3d sumVelInfl = Eigen::Vector3d::Zero();
        
        if (accelerate && timeStep >0){
            FMM.barnesHut((*bPanels)[i]->getCenter());
        }
        
        for(int j=0; j<filaments.size(); j++)
        {
            sumVelInfl+=filaments[j]->velInfl((*bPanels)[i]->getCenter());
        }
        (*bPanels)[i]->setSigma((Vinf+sumVelInfl),0);
        sigmas(i) = (*bPanels)[i]->getSigma();
    }
}

bool cpCase::solveMatrixEq()
{
    bool converged = true;
    
    Eigen::MatrixXd* A = geom->getA();
    Eigen::MatrixXd* B = geom->getB();
    Eigen::VectorXd RHS = -(*B)*sigmas;
    Eigen::VectorXd doubletStrengths(bPanels->size());

    
    Eigen::BiCGSTAB<Eigen::MatrixXd> res;
    res.compute((*A));
    doubletStrengths = res.solve(RHS);
    if (res.error() > pow(10,-10))
    {
        converged = false;
    }
    
    for (int i=0; i<bPanels->size(); i++)
    {
        (*bPanels)[i]->setMu(doubletStrengths(i));
        (*bPanels)[i]->setPotential(Vinf);
    }
    
    for (int i=0; i<wPanels->size(); i++)
    {
        (*wPanels)[i]->setMu();
        (*wPanels)[i]->setPotential(Vinf);
    }
    
    return converged;
}

bool cpCase::solveVPmatrixEq()
{
    bool converged = true;
    
    Eigen::MatrixXd* A = geom->getA();
    Eigen::MatrixXd* B = geom->getB();
    Eigen::VectorXd RHS = -(*B)*(sigmas);
    Eigen::VectorXd doubletStrengths(bPanels->size());
    
    
    Eigen::BiCGSTAB<Eigen::MatrixXd> res;
    res.compute((*A));
    doubletStrengths = res.solve(RHS);
    if (res.error() > pow(10,-10))
    {
        converged = false;
    }
    
    for (int i=0; i<bPanels->size(); i++)
    {
        (*bPanels)[i]->setMu(doubletStrengths(i));
        (*bPanels)[i]->setPotential(Vinf);

    }
    for (int i=0; i<wPanels->size(); i++)
    {
        (*wPanels)[i]->setPrevStrength((*wPanels)[i]->getMu());
        (*wPanels)[i]->setMu();
        (*wPanels)[i]->setPotential(Vinf);
    }
    return converged;
}


void cpCase::compVelocity()
{
    //  Velocity Survey with known doublet and source strengths
    CM.setZero();
    Eigen::Vector3d moment;
    Fbody = Eigen::Vector3d::Zero();
    bodyPanel* p;
    
    
    for (int i=0; i<bPanels->size(); i++)
    {
        p = (*bPanels)[i];
        Eigen::Vector3d sumPartInfl = Eigen::Vector3d::Zero();
        for(int j=0; j<particles.size(); j++){
            sumPartInfl += particles[j]->partVelInfl(p->getCenter());
        }
        p->computeVelocity(PG,(Vinf),sumPartInfl);
        p->computeCp(Vmag);
        Fbody += -p->getCp()*p->getArea()*p->getBezNormal()/params->Sref;
        moment = p->computeMoments(params->cg);
        CM(0) += moment(0)/(params->Sref*params->bref);
        CM(1) += moment(1)/(params->Sref*params->cref);
        CM(2) += moment(2)/(params->Sref*params->bref);
    }
    Fwind = bodyToWind(Fbody);
//    std::cout << "CL = " << Fbody.z() << std::endl;
    CL.push_back(Fbody.z());
    
}

void cpCase::trefftzPlaneAnalysis()
{
    std::vector<wake*> wakes = geom->getWakes();
    CL_trefftz = 0;
    CD_trefftz = 0;
    double CD_trefftzVel = 0;
    for (int i=0; i<wakes.size(); i++)
    {
        wakes[i]->trefftzPlane(Vmag,params->Sref);
        CD_trefftzVel += wakes[i]->trefftzPlaneFromVel(Vmag, params->Sref);
        CL_trefftz += wakes[i]->getCL()/PG;
        CD_trefftz += wakes[i]->getCD()/pow(PG,2);
    }
    std::cout << "CD from Velocity trefftz plane: " << CD_trefftzVel << std::endl;
}

void cpCase::createStreamlines()
{
    // Gather TE Panels
    std::vector<surface*> surfs = geom->getSurfaces();
    std::vector<std::pair<Eigen::Vector3d,bodyPanel*>> streamPnts;
    bodyStreamline* s;
    
    for (int i=0; i<surfs.size(); i++)
    {
        streamPnts = surfs[i]->getStreamlineStartPnts(Vinf,PG);
        for (int j=0; j<streamPnts.size(); j++)
        {
            s = new bodyStreamline(std::get<0>(streamPnts[j]),std::get<1>(streamPnts[j]),Vinf,PG,geom,3,false);
            bStreamlines.push_back(s);
        }
    }
    
//    // Off Body Streamline (Temporary)
//    std::vector<streamline*> sLines;
//    streamline* sline;
//    Eigen::Vector3d start;
//    int ny = 1;
//    int nz = 6;
//    double dz = 0.1;
//    for (int i=0; i<ny; i++)
//    {
//        for (int j=0; j<nz; j++)
//        {
////            start << -1,-0.5*params->bref+params->bref*i/(ny-1),-0.35+j*dz;
//            start << -1,-2.85,-0.45+j*dz;
//            sline = new streamline(start,20,0.0001,Vinf,PG,geom);
//            sLines.push_back(sline);
//        }
//    }
//    
//    piece p;
//    std::vector<piece> pieces;
//    std::vector<pntDataArray> data;
//    pntDataArray vel("Velocity");
//    Eigen::MatrixXi con;
//    Eigen::MatrixXd pntMat;
//    std::vector<Eigen::Vector3d> pnts,velocities;
//    
//    for (int i=0; i<sLines.size(); i++)
//    {
//        pnts = sLines[i]->getPnts();
//        velocities = sLines[i]->getVelocities();
//        vel.data.resize(velocities.size(),3);
//        pntMat.resize(pnts.size(),3);
//        con.resize(pnts.size()-1,2);
//        for (int j=0; j<pnts.size(); j++)
//        {
//            pntMat.row(j) = pnts[j];
//            vel.data.row(j) = velocities[j];
//            if (j<con.rows())
//            {
//                con(j,0) = j;
//                con(j,1) = j+1;
//            }
//        }
//        data.push_back(vel);
//        p.pnts = pntMat;
//        p.connectivity = con;
//        p.pntData = data;
//        
//        pieces.push_back(p);
//        data.clear();
//    }
//    
//    boost::filesystem::path path = boost::filesystem::current_path();
//    
//    std::string fname = path.string()+"/OffBodyStreamlines.vtu";
//    VTUfile wakeFile(fname,pieces);
//    
//    for (int i=0; i<sLines.size(); i++)
//    {
//        delete sLines[i];
//    }
}

void cpCase::stabilityDerivatives()
{
    double delta = 0.5;
    double dRad = delta*M_PI/180;
    cpCase dA(geom,Vmag,alpha+delta,beta,mach,params);
    cpCase dB(geom,Vmag,alpha,beta+delta,mach,params);
    
    dA.run(false,false,false,false);
    dB.run(false,false,false,false);
    
    Eigen::Vector3d FA = dA.getWindForces();
    FA(2) = dA.getCL();
    FA(0) = dA.getCD();
    
    Eigen::Vector3d FB = dB.getWindForces();
    FB(2) = dB.getCL();
    FB(0) = dB.getCD();
    
    Eigen::Vector3d F = Fwind;
    F(2) = CL_trefftz;
    F(0) = CD_trefftz;
    
    // Finish calculating stability derivatives
    dF_dAlpha = (FA-F)/dRad;
    dF_dBeta = (FB-F)/dRad;
    
    dM_dAlpha = (dA.getMoment()-CM)/dRad;
    dM_dBeta = (dB.getMoment()-CM)/dRad;
    
    
}

void cpCase::writeFiles()
{
    std::stringstream caseLabel;
    caseLabel << "/V" << Vmag << "_Mach" << mach << "_alpha" << alpha << "_beta" << beta;
    boost::filesystem::path subdir = boost::filesystem::current_path().string()+caseLabel.str();
    if (!boost::filesystem::exists(subdir))
    {
        boost::filesystem::create_directories(subdir);
    }
    Eigen::MatrixXd nodeMat = geom->getNodePnts();
    writeBodyData(subdir,nodeMat);
    if (geom->getWakes().size() > 0)
    {
        writeWakeData(subdir,nodeMat);
        writeSpanwiseData(subdir);
    }
    
    
    if(vortPartFlag){
        if(timeStep > 0){
            writeParticleData(subdir);
            writeFilamentData(subdir);
        }
    }
    
    if (params->surfStreamFlag)
    {
        writeBodyStreamlines(subdir);
    }

}

void cpCase::writeBodyData(boost::filesystem::path path,const Eigen::MatrixXd &nodeMat)
{
    std::vector<cellDataArray> data;
    cellDataArray mu("Doublet Strengths"),sigma("Source Strengths"),pot("Velocity Potential"),V("Velocity"),Cp("Cp"),bN("bezNormals"),x("xPosition"),y("yPosition"),z("zPostition");
    Eigen::MatrixXi con(bPanels->size(),3);
    mu.data.resize(bPanels->size(),1);
    sigma.data.resize(bPanels->size(),1);
    pot.data.resize(bPanels->size(),1);
    V.data.resize(bPanels->size(),3);
    Cp.data.resize(bPanels->size(),1);
    bN.data.resize(bPanels->size(),3);
    x.data.resize(bPanels->size(),1);
    y.data.resize(bPanels->size(),1);
    z.data.resize(bPanels->size(),1);
    
    for (int i=0; i<bPanels->size(); i++)
    {
        mu.data(i,0) = (*bPanels)[i]->getMu();
        sigma.data(i,0) = (*bPanels)[i]->getSigma();
        pot.data(i,0) = (*bPanels)[i]->getPotential();
        V.data.row(i) = (*bPanels)[i]->getGlobalV();
        Cp.data(i,0) = (*bPanels)[i]->getCp();
        con.row(i) = (*bPanels)[i]->getVerts();
        bN.data.row(i) = (*bPanels)[i]->getBezNormal();
        x.data(i,0) = (*bPanels)[i]->getCenter().x();
        y.data(i,0) = (*bPanels)[i]->getCenter().y();
        z.data(i,0) = (*bPanels)[i]->getCenter().z();
    }
    
    data.push_back(mu);
    data.push_back(sigma);
    data.push_back(pot);
    data.push_back(V);
    data.push_back(Cp);
    data.push_back(bN);
    data.push_back(x);
    data.push_back(y);
    data.push_back(z);
    
    piece body;
    body.pnts = nodeMat;
    body.connectivity = con;
    body.cellData = data;
    
    std::string fname = path.string()+"/surfaceData-" + std::to_string(timeStep)+".vtu";
    VTUfile bodyFile(fname,body);
}

void cpCase::writeWakeData(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat)
{
    std::vector<cellDataArray> data;
    cellDataArray mu("Doublet Strengths"),pot("Velocity Potential");
    Eigen::MatrixXi con;
    if(vortPartFlag){ //VPP
        con.resize(wPanels->size(),4);
    }else{
        con.resize(wPanels->size(),3);
    }
    mu.data.resize(wPanels->size(),1);
    pot.data.resize(wPanels->size(),1);
    for (int i=0; i<wPanels->size(); i++)
    {
        mu.data(i,0) = (*wPanels)[i]->getMu();
        pot.data(i,0) = (*wPanels)[i]->getPotential();
        con.row(i) = (*wPanels)[i]->getVerts();
    }
    data.push_back(mu);
    data.push_back(pot);
    
    piece wake;
    wake.pnts = nodeMat;
    wake.connectivity = con;
    wake.cellData = data;
    std::string fname = path.string()+"/wakeData-"+std::to_string(timeStep)+".vtu";
    VTUfile wakeFile(fname,wake);
}

void cpCase::writeFilamentData(boost::filesystem::path path)
{
    std::vector<cellDataArray> data;
    cellDataArray mu("Gamma");
    Eigen::MatrixXi con;
    
    Eigen::MatrixXd nodeMat(2*filaments.size(),3);
    for(int i=0; i<filaments.size(); i++){
        nodeMat.row(2*i) = filaments[i]->getP1();
        nodeMat.row(2*i+1) = filaments[i]->getP2();
    }
    

    con.resize(filaments.size(),2); // 2 pnts per filament
    
    mu.data.resize(filaments.size(),1);

    int nodeCounter = 0; // Used to make filament end points
    for (int i=0; i<filaments.size(); i++)
    {
        mu.data(i,0) = filaments[i]->getStrength();
        con(i,0) = nodeCounter;
        con(i,1) = nodeCounter+1;
        nodeCounter = nodeCounter+2;
    }
    
    data.push_back(mu);
    
    piece fils;
    fils.pnts = nodeMat;
    fils.connectivity = con;
    fils.cellData = data;
    std::string fname = path.string()+"/filaments-"+std::to_string(timeStep)+".vtu";
    VTUfile filFile(fname,fils);
}

void cpCase::writeParticleData(boost::filesystem::path path)
{

    Eigen::MatrixXd partMat(particles.size(),3);
    for (int i=0; i<particles.size(); i++)
    {
        partMat.row(i) = particles[i]->pos;
    }
    
    std::vector<cellDataArray> data;
    cellDataArray strength("Strength");
    Eigen::MatrixXi con(particles.size(),1);
    
    strength.data.resize(particles.size(),3);
    for (int i=0; i<particles.size(); i++)
    {
        strength.data.row(i) = particles[i]->strength;
        con(i) = i;
    }
    data.push_back(strength);
    
    piece parts;
    parts.pnts = partMat;
    parts.connectivity = con;
    parts.cellData = data;
    
    std::string fname = path.string()+"/particleData-"+std::to_string(timeStep)+".vtu";
    VTUfile partFile(fname,parts);
}


void cpCase::writeSpanwiseData(boost::filesystem::path path)
{
    std::vector<wake*> wakes = geom->getWakes();
    for (int i=0; i<wakes.size(); i++)
    {
        Eigen::VectorXd spanLoc,Cl,Cd;
        spanLoc = 2*wakes[i]->getSpanwisePnts()/params->bref;
        Cl = wakes[i]->getSpanwiseCl()/PG;
        Cd = wakes[i]->getSpanwiseCd()/pow(PG,2);
        
        std::stringstream ss;
        ss << path.string() << "/spanwiseData_Wake" << i+1 << ".csv";
        std::string fname = ss.str();
        std::ofstream fout;
        fout.open(fname);
        if (fout)
        {
            fout << "2y/b,Cl,Cdi" << std::endl;
            for (int i=0; i<spanLoc.size(); i++)
            {
                fout << spanLoc(i) << "," << Cl(i) << "," << Cd(i) << std::endl;
            }
        }
        fout.close();
    }
}

void cpCase::writeBodyStreamlines(boost::filesystem::path path)
{
    piece p;
    std::vector<piece> pieces;
    std::vector<pntDataArray> data;
    pntDataArray vel("Velocity");
    Eigen::MatrixXi con;
    Eigen::MatrixXd pntMat;
    std::vector<Eigen::Vector3d> pnts,velocities;
    
    for (int i=0; i<bStreamlines.size(); i++)
    {
        pnts = bStreamlines[i]->getPnts();
        velocities = bStreamlines[i]->getVelocities();
        vel.data.resize(velocities.size(),3);
        pntMat.resize(pnts.size(),3);
        con.resize(pnts.size()-1,2);
        for (int j=0; j<pnts.size(); j++)
        {
            pntMat.row(j) = pnts[j];
            vel.data.row(j) = velocities[j];
            if (j<con.rows())
            {
                con(j,0) = j;
                con(j,1) = j+1;
            }
        }
        data.push_back(vel);
        p.pnts = pntMat;
        p.connectivity = con;
        p.pntData = data;
        
        pieces.push_back(p);
        data.clear();
    }
    
    std::string fname = path.string()+"/streamlines.vtu";
    VTUfile wakeFile(fname,pieces);
}

void cpCase::collapseWakeForEachEdge(){
    // Go through each panel. Ignore the trialing edge as it has no circulation. After creating a particle at each edge, put a pointer to the edge in a vector. Before making any particles, make sure the edge has not yet been used. Edge collapse follows
    
    std::vector<edge*> usedEdges;
    for(int i=0; i<(*wPanels).size(); i++){
        std::vector<edge*> pEdges = (*wPanels)[i]->edgesInOrder();
        for (int j=1; j<4; j++) {
            if (!edgeIsUsed(pEdges[j],usedEdges)) {
                usedEdges.push_back(pEdges[j]);
                Eigen::Vector3d pos = seedPos((*wPanels)[i], j);
                Eigen::Vector3d strength = edgeStrength((*wPanels)[i], pEdges[j], j);
                double radius = (*wPanels)[i]->getPartRadius(Vinf,dt);
                
                particle* p = new particle(pos, strength, radius, {0,0,0}, {0,0,0}); //last two are previous pos and strength values used for advanced time stepper.
                particles.push_back(p);
            }
        }
    }

    // Create filament
    if(timeStep == 1){
        for(int i=0; i<(*wPanels).size(); i++){
            vortexFil* fil;
            Eigen::Vector3d p1,p2;
            p1 = (*wPanels)[i]->pointsInOrder()[2]->getPnt();
            p2 = (*wPanels)[i]->pointsInOrder()[3]->getPnt();
            
            if(startingWake){
                fil = new vortexFil(p1, p2,-(*wPanels)[i]->getMu(), (*wPanels)[i]); // Negative strength is because filament is actually the upstream edge being convected which is oriented the opposite direction as downstream edge
            }else{
                fil = new vortexFil(p1, p2,0, (*wPanels)[i]); // Negative strength is because filament is actually the upstream edge being convected which is oriented the opposite direction as downstream edge
            }
            filaments.push_back(fil);
            (*wPanels)[i]->setVortFil(fil);
        }
    }else{
        if(startingWake){
            for(int i=0; i<(*wPanels).size(); i++){
                filaments[i]->setStrength(-(*wPanels)[i]->getMu()); // Negative strength is because filament is actually the upstream edge being convected which is oriented the opposite direction as downstream edge
            }
        }else{
            for(int i=0; i<(*wPanels).size(); i++){
                filaments[i]->setStrength(-(*wPanels)[i]->getMu());
            }
        }
    }
    
}

bool cpCase::edgeIsUsed(edge* thisEdge, std::vector<edge*> pEdges){
    
    for(int i=0; i<pEdges.size(); i++){
        if(thisEdge == pEdges[i]){
            return true;
        }
    }
    return false;
}

Eigen::Vector3d cpCase::edgeStrength(wakePanel* pan, edge* curEdge, int edgeNum){

    Eigen::Vector3d strength;
    
//    if(edgeNum == 0){
//        // Equation for part strength comes from Martin eq. 5.25
//        Eigen::Vector3d Rj = pan->pointsInOrder()[0]->getPnt();
//        Eigen::Vector3d Ri = pan->pointsInOrder()[1]->getPnt();
//        strength = (pan->getMu())*(Ri-Rj);
    //    }
    if(edgeNum == 2) { // If is far edge, don't need to worry about neighbor panel
        Eigen::Vector3d Rj = pan->pointsInOrder()[2]->getPnt();
        Eigen::Vector3d Ri = pan->pointsInOrder()[3]->getPnt();
        if (startingWake) {
            strength = (pan->getMu()-pan->getPrevStrength())*(Ri-Rj);
        }
        else{
//            if(timeStep == 1){
//                strength = Eigen::Vector3d::Zero();
//            }else if(timeStep == 2){
//                strength = (pan->getMu())*(Ri-Rj);
//            }else{
//                strength = (pan->getMu()-pan->getPrevStrength())*(Ri-Rj);
//            }
            strength = Eigen::Vector3d::Zero();

        }
    }
    else if(edgeNum == 1){
        wakePanel* otherPan = curEdge->getOtherWakePan(pan);
        Eigen::Vector3d Rj = pan->pointsInOrder()[1]->getPnt();
        Eigen::Vector3d Ri = pan->pointsInOrder()[2]->getPnt();
        if(otherPan){ // Panel has neighbor
            strength = (pan->getMu()-otherPan->getMu())*(Ri-Rj);
        }else
            strength = pan->getMu()*(Ri-Rj);
    }
    else{
        wakePanel* otherPan = curEdge->getOtherWakePan(pan);
        Eigen::Vector3d Rj = pan->pointsInOrder()[3]->getPnt();
        Eigen::Vector3d Ri = pan->pointsInOrder()[0]->getPnt();
        if(otherPan){ // Panel has neighbor
            strength = (pan->getMu()-otherPan->getMu())*(Ri-Rj);
        }else
            strength = pan->getMu()*(Ri-Rj);
    }
    return strength;
}

Eigen::Vector3d cpCase::seedPos(wakePanel* pan, int edgeNum){
    
    std::vector<cpNode*> nodes = pan->pointsInOrder();
    
    if(edgeNum == 0){      // Avg of first proj nodes for pntsIO[0] and pntsIO[1]
        return (nodes[0]->firstProjNode(dt, Vinf.norm()) + nodes[1]->firstProjNode(dt, Vinf.norm()))/2;
    }
    else if(edgeNum == 1){ // Avg of 1st and second proj nodes of pntsIO[1]
        return (nodes[1]->firstProjNode(dt, Vinf.norm()) + nodes[1]->secProjNode(dt, Vinf.norm()))/2;
    }
    else if(edgeNum == 2){ // Avg of second proj nodes for pntsIO[0] and pntsIO[1]
        return (nodes[0]->secProjNode(dt, Vinf.norm()) + nodes[1]->secProjNode(dt, Vinf.norm()))/2;
    }
    else if(edgeNum == 3){ // Avg of 1st and second proj nodes of pntsIO[0]
        return (nodes[0]->firstProjNode(dt, Vinf.norm()) + nodes[0]->secProjNode(dt, Vinf.norm()))/2;
    }
    
    std::cout << "Error: wrong edge number for particle seed position..." << std::endl;
    std::exit(0);
    
}

Eigen::Vector3d cpCase::velocityInflFromEverything(Eigen::Vector3d POI){
        
    // Freestream influence
    Eigen::Vector3d velOnPart = Vinf;
    
    // Body panel influence
    for(int j=0;j<(*bPanels).size();j++){
        velOnPart += (*bPanels)[j]->panelV(POI);
    }
    
    // Buffer wake influence
    for(int j=0;j<(*wPanels).size();j++){
        velOnPart += (*wPanels)[j]->panelV(POI);
    }
    
    
    // Particle influence
    if(accelerate){
        velOnPart += FMM.barnesHut(POI);
    }
    else{
        for(int j=0;j<particles.size();j++)
        {
            velOnPart += particles[j]->partVelInfl(POI);
        }
    }

    // Vortex Filament influence
    for(int i=0; i<filaments.size(); i++){
        velOnPart +=filaments[i]->velInfl(POI);
    }
    
    return velOnPart;
}

void cpCase::convectParticles(){
    std::vector<Eigen::Vector3d> newPartPositions;
    
//    particleFMM pFMM;
//    pFMM.build(&partOctree);

    
    for(int i=0;i<particles.size();i++){
        
        Eigen::Vector3d newPos;
        
        if(rungeKuttaStep){
            Eigen::Vector3d POI = particles[i]->pos;
            Eigen::Vector3d k1 = velocityInflFromEverything(POI);
            Eigen::Vector3d k2 = velocityInflFromEverything(POI+k1*dt/2);
            Eigen::Vector3d k3 = velocityInflFromEverything(POI+k2*dt/2);
            Eigen::Vector3d k4 = velocityInflFromEverything(POI+k3*dt);

            newPartPositions.push_back(POI + dt*(k1/6 + k2/3 + k3/3 + k4/6));
            
        }
        else // Do Adams Bashforth
        {
            Eigen::Vector3d velOnPart = velocityInflFromEverything(particles[i]->pos);
            
            if(particles[i]->getPrevVelInfl().isZero())
            {
                newPos = particles[i]->pos + dt*velOnPart;
            }
            else{
                // Adams bashforth scheme
                newPos = particles[i]->pos + dt*(1.5*velOnPart - 0.5*particles[i]->getPrevVelInfl());
            }

            particles[i]->setPrevVelInfl(velOnPart);
            newPartPositions.push_back(newPos);
        }
    }
    
    for(int i=0;i<particles.size();i++){
        particles[i]->setPos(newPartPositions[i]);

    }
}



void cpCase::particleStrengthUpdate(){
    // This function uses the combined vortex stretching and diffusion equation used by Wincklemans for a regularized vortex core with high algebraic smoothing. (he refers to it as the strength update equation.)
    
    std::vector<Eigen::Vector3d> stretchDiffVec; // Creating a vector of diffusion values because the strength change needs to be set after all particle influences have been calculated
    for(int i=0; i<particles.size(); i++)
    {
        Eigen::Vector3d dAlpha = Eigen::Vector3d::Zero();
        
        // Stretching from particles
        for(int j=0; j<particles.size(); j++)
        {
                dAlpha += particles[i]->partStrengthUpdate(particles[j]);
        }
        
        // Stretching from body panels
        for(int j=0; j<(*bPanels).size(); j++)
        {
            dAlpha += (*bPanels)[j]->partStretching(particles[i]);
        }
        
        // Stretching from wake panels
        for(int j=0;j<(*wPanels).size(); j++)
        {
            dAlpha += (*wPanels)[j]->partStretching(particles[i]);
        }
        
        stretchDiffVec.push_back(dAlpha);
    }
    
    // No need for Kutta accuracy
    for(int i=0;i<particles.size();i++){
        Eigen::Vector3d newStrength;
        if(particles[i]->getprevStrengthUpdate().isZero())
        {
            newStrength = particles[i]->strength + stretchDiffVec[i]*dt;
        }else{
            //Adams bashforth
            newStrength = particles[i]->strength + dt*(1.5*stretchDiffVec[i] - 0.5*particles[i]->getprevStrengthUpdate());
        }
        particles[i]->setprevStrengthUpdate(stretchDiffVec[i]);
        particles[i]->setStrength(newStrength);
    }
    
}

void cpCase::particleStrengthUpdateGaussian(){
    // This function uses the update equations found in 'Vortex Methods for DNS of...' by Plouhmhans. It uses the particle strength exchange for the viscous diffusion
    
    std::vector<Eigen::Vector3d> stretchDiffVec; // Creating vector values because the strength change needs to be set after all particle influences have been calculated
    for(int i=0; i<particles.size(); i++){
        
        Eigen::Vector3d dAlpha_diff = Eigen::Vector3d::Zero();
        Eigen::Vector3d dAlpha_stretch = Eigen::Vector3d::Zero();
        for(int j=0; j<particles.size(); j++){
            if(i!=j){  // Shouldnt need this here! test and take out!
                dAlpha_diff += particles[i]->viscousDiffusionGaussian(particles[j]);
                dAlpha_stretch += particles[i]->vortexStretchingGaussian(particles[j]);
            }
        }
        stretchDiffVec.push_back(dAlpha_diff+dAlpha_stretch);
    }
    
    for(int i=0;i<particles.size();i++){
        Eigen::Vector3d newStrength;
        if(particles[i]->getprevStrengthUpdate().isZero())
        {
            newStrength = particles[i]->strength + stretchDiffVec[i]*dt;
        }else
        {
            newStrength = particles[i]->strength + dt*(1.5*stretchDiffVec[i] - 0.5*particles[i]->getprevStrengthUpdate()); //Adams bashforth
        }
        particles[i]->setprevStrengthUpdate(stretchDiffVec[i]);
        particles[i]->setStrength(newStrength);
    }
}

//double cpCase::trefftzPlaneCd(std::vector<particle*> particles){
//    int yPnts = 30;
//    int zPnts = 30;
//    double xPartMax=particles[0]->pos.x();
//    double yPartMax=particles[0]->pos.y(); // initializing with first particle so that
//    double yPartMin=particles[0]->pos.y();
//    double zPartMax=particles[0]->pos.z();
//    double zPartMin=particles[0]->pos.z();
//    
//    for(int i=0; i<particles.size(); i++){
//        if(particles[i]->pos.x() > xPartMax){
//            xPartMax = particles[i]->pos.x();
//        }
//        if(particles[i]->pos.y() > yPartMax){
//            yPartMax = particles[i]->pos.y();
//        }
//        if(particles[i]->pos.y() < yPartMin){
//            yPartMin = particles[i]->pos.y();
//        }
//        if(particles[i]->pos.z() > zPartMax){
//            zPartMax = particles[i]->pos.z();
//        }
//        if(particles[i]->pos.z() < zPartMax){
//            zPartMin = particles[i]->pos.z();
//        }
//    }
//
//    Eigen::MatrixXd w,v,Cd;
//    Eigen::Vector3d pnt;
//    double yLoc,zLoc;
//    w.resize(zPnts,yPnts);
//    v.resize(zPnts,yPnts);
//    Cd.resize(zPnts,yPnts);
//    double yStep = (yPartMax-yPartMin)/(yPnts);
//    double zStep = (zPartMax-zPartMin)/(zPnts);
//
//    std::cout << "Calculating particle Trefftz Plane" << std::endl;
//    for(int i=0; i<yPnts; i++){
//        yLoc = yPartMin+i*yStep;
//        for(int j=0; j<zPnts; j++){
//            zLoc = zPartMin+zStep;
//            pnt << xPartMax,yLoc,zLoc;
//            for(int k=0; k<particles.size(); k++){
//                v(i,j) += particles[k]->partVelInfl(pnt).y();
//                w(i,j) += particles[k]->partVelInfl(pnt).z();
//            }
//            Cd(i,j) = (v(i,j)*v(i,j)+w(i,j)*w(i,j));
//        }
//    }
//    
//    return Cd.sum();
//
//
////    for (int i=1; i<nPnts; i++)
////    {
////        yLoc(i) = yMin+i*step;
////        pWake = pntInWake(xTrefftz, yLoc(i));
////        Cl(i) = 2*dPhi(i)/(Vinf*Sref);
////        Cd(i) = dPhi(i)*w(i)/(Vinf*Vinf*Sref);
////    }
//};
//
//


//void cpCase::collapseBufferWake(){
//    std::cout << "collapse buffer wake needs to be modified to include correct edge strengths and vortex filaments" << std::endl;
//
//    for(int i=0;i<(*wPanels).size();i++){
//
//        Eigen::Vector3d pos = (*wPanels)[i]->partSeedPt(Vinf,dt);
//        double radius = (*wPanels)[i]->getPartRadius(Vinf,dt);
//        Eigen::Vector3d strength;
//        if(timeStep == 1){
//            strength = (*wPanels)[i]->panToPartStrengthT1();
////            strength = (*wPanels)[i]->panToPartStrength();
//        }else{
//            strength = (*wPanels)[i]->panToPartStrength();
//        }
//
//        particle* p = new particle(pos, strength, radius, {0,0,0}, {0,0,0}); // last two are previous positon and strength values used for advanced time stepper
//        particles.push_back(p);
//    }
//}
