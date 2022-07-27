// Copyright (c) 2009-2019 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.


// Maintainer: joaander

#include "TwoStepBDGPU.h"
#include "TwoStepBDGPU.cuh"

#ifdef ENABLE_MPI
#include "hoomd/HOOMDMPI.h"
#endif

#include <fstream>
#include <iostream>

namespace py = pybind11;

using namespace std;

/*! \file TwoStepBDGPU.h
    \brief Contains code for the TwoStepBDGPU class
*/

/*! \param sysdef SystemDefinition this method will act on. Must not be NULL.
    \param group The group of particles this integration method is to work on
    \param T Temperature set point as a function of time
    \param seed Random seed to use in generating random numbers
    \param use_lambda If true, gamma=lambda*diameter, otherwise use a per-type gamma via setGamma()
    \param lambda Scale factor to convert diameter to gamma
*/
TwoStepBDGPU::TwoStepBDGPU(std::shared_ptr<SystemDefinition> sysdef,
                           std::shared_ptr<ParticleGroup> group,
                           std::shared_ptr<Variant> T,
                           unsigned int seed,
                           bool use_lambda,
                           Scalar lambda,
                           bool noiseless_t,
                           bool noiseless_r)
    : TwoStepBD(sysdef, group, T, seed, use_lambda, lambda, noiseless_t, noiseless_r)
    {

    pitsFile.open(name);
    x_pitCPU = new double[SIZE];                    //vidilenie pam9ti dl9 cpu massiva
    y_pitCPU = new double[SIZE];    
    
    cudaMalloc(&x_pit, SIZE*sizeof(double));        //vididlenie pam9ti dl9 gpu massiva
    cudaMalloc(&y_pit, SIZE*sizeof(double));
    
    //4etaem file s koordinatami potentialnih 9m
    pitsFile >> epsilon1;                                                        //glubina
    pitsFile >> mySigma1;                                                        //radius
    mySigma1 *= mySigma1;
        
    while(nMax < SIZE){
        
        double temp1 = 0;
        double temp2 = 0;
            
        pitsFile >> temp1;
        pitsFile >> temp2;
            
        if(temp1 == 1337228 && temp2 == 1488322)
            break;
            
        else{
                
            x_pitCPU[nMax1] = temp1;
            y_pitCPU[nMax1] = temp2;
            nMax1++;
                
        }
            
    }
    pitsFile.close();
    
    for(int i = 0 ; i < nMax1 ; i++)
        std::cout << x_pitCPU[i] << "\t" << y_pitCPU[i] << std::endl;
    
    //PERENOS DANNIH S CPU NA GPU
    cudaMemcpy(x_pit, x_pitCPU, nMax1*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(y_pit, y_pitCPU, nMax1*sizeof(double), cudaMemcpyHostToDevice);
    
    if (!m_exec_conf->isCUDAEnabled())
        {
        m_exec_conf->msg->error() << "Creating a TwoStepBDGPU while CUDA is disabled" << endl;
        throw std::runtime_error("Error initializing TwoStepBDGPU");
        }

    m_block_size = 256;
    }

/*! \param timestep Current time step
    \post Particle positions are moved forward a full time step and velocities are redrawn from the proper distribution.
*/

TwoStepBDGPU::~TwoStepBDGPU(){
    
    delete [] x_pitCPU;
    delete [] y_pitCPU;
    
    cudaFree(x_pit);
    cudaFree(y_pit);
    
}

void TwoStepBDGPU::integrateStepOne(unsigned int timestep)
    {
    
    // profile this step
    if (m_prof)
        m_prof->push(m_exec_conf, "BD step 1");

    // access all the needed data
    BoxDim box = m_pdata->getBox();
    ArrayHandle< unsigned int > d_index_array(m_group->getIndexArray(), access_location::device, access_mode::read);
    unsigned int group_size = m_group->getNumMembers();
    const unsigned int D = m_sysdef->getNDimensions();
    const GlobalArray< Scalar4 >& net_force = m_pdata->getNetForce();

    ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar4> d_vel(m_pdata->getVelocities(), access_location::device, access_mode::readwrite);
    ArrayHandle<int3> d_image(m_pdata->getImages(), access_location::device, access_mode::readwrite);
    
    //dobavil massiv s zna4eni9mi potntialnih 9m
    //ArrayHandle<Scalar4> r_pits(m_pdata->getVelocities(), access_location::device, access_mode::readwrite);

    ArrayHandle<Scalar4> d_net_force(net_force, access_location::device, access_mode::read);
    ArrayHandle<Scalar> d_gamma(m_gamma, access_location::device, access_mode::read);
    ArrayHandle<Scalar> d_diameter(m_pdata->getDiameters(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_tag(m_pdata->getTags(), access_location::device, access_mode::read);

    // for rotational noise
    ArrayHandle<Scalar3> d_gamma_r(m_gamma_r, access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_orientation(m_pdata->getOrientationArray(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar4> d_torque(m_pdata->getNetTorqueArray(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar3> d_inertia(m_pdata->getMomentsOfInertiaArray(), access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_angmom(m_pdata->getAngularMomentumArray(), access_location::device, access_mode::readwrite);
    
    langevin_step_two_args args;
    args.d_gamma = d_gamma.data;
    args.n_types = m_gamma.getNumElements();
    args.use_lambda = m_use_lambda;
    args.lambda = m_lambda;
    args.T = m_T->getValue(timestep);
    args.timestep = timestep;
    args.seed = m_seed;
    args.d_sum_bdenergy = NULL;
    args.d_partial_sum_bdenergy = NULL;
    args.block_size = m_block_size;
    args.num_blocks = 0; // handled in driver function
    args.tally = false;

    bool aniso = m_aniso;

    if (m_exec_conf->allConcurrentManagedAccess())
        {
        // prefetch gammas
        auto& gpu_map = m_exec_conf->getGPUIds();
        for (unsigned int idev = 0; idev < m_exec_conf->getNumActiveGPUs(); ++idev)
            {
            cudaMemPrefetchAsync(m_gamma.get(), sizeof(Scalar)*m_gamma.getNumElements(), gpu_map[idev]);
            cudaMemPrefetchAsync(m_gamma_r.get(), sizeof(Scalar)*m_gamma_r.getNumElements(), gpu_map[idev]);
            }
        if (m_exec_conf->isCUDAErrorCheckingEnabled())
            CHECK_CUDA_ERROR();
        }

    m_exec_conf->beginMultiGPU();
    
    // perform the update on the GPU
    gpu_brownian_step_one(d_pos.data,                               
                          d_vel.data,
                          x_pit,                                        //PEREDAU NUJNIE ZNACHENI9
                          y_pit,                                        //
                          mySigma1,                                       //
                          epsilon1,                                       //
                          nMax1,                                          //
                          d_image.data,
                          box,
                          d_diameter.data,
                          d_tag.data,
                          d_index_array.data,
                          group_size,
                          d_net_force.data,
                          d_gamma_r.data,
                          d_orientation.data,
                          d_torque.data,
                          d_inertia.data,
                          d_angmom.data,
                          args,
                          aniso,
                          m_deltaT,
                          D,
                          m_noiseless_t,
                          m_noiseless_r,
                          m_group->getGPUPartition());

    if(m_exec_conf->isCUDAErrorCheckingEnabled())
        CHECK_CUDA_ERROR();

    m_exec_conf->endMultiGPU();

    // done profiling
    if (m_prof)
        m_prof->pop(m_exec_conf);
    
    }

/*! \param timestep Current time step
    \post particle velocities are moved forward to timestep+1 on the GPU
*/
void TwoStepBDGPU::integrateStepTwo(unsigned int timestep)
    {
    // there is no step 2
    }

void export_TwoStepBDGPU(py::module& m)
    {
    py::class_<TwoStepBDGPU, std::shared_ptr<TwoStepBDGPU> >(m, "TwoStepBDGPU", py::base<TwoStepBD>())
        .def(py::init< std::shared_ptr<SystemDefinition>,
                               std::shared_ptr<ParticleGroup>,
                               std::shared_ptr<Variant>,
                               unsigned int,
                               bool,
                               Scalar,
                               bool,
                               bool>())
        ;
    }
