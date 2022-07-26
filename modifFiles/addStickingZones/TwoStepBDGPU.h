// Copyright (c) 2009-2019 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.


// Maintainer: joaander

#include "TwoStepBD.h"

#ifndef __TWO_STEP_BD_GPU_H__
#define __TWO_STEP_BD_GPU_H__

/*! \file TwoStepBDGPU.h
    \brief Declares the TwoStepBDGPU class
*/

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

#include <hoomd/extern/pybind/include/pybind11/pybind11.h>
#include <stdio.h>

//! Implements Brownian dynamics on the GPU
/*! GPU accelerated version of TwoStepBD

    \ingroup updaters
*/
class PYBIND11_EXPORT TwoStepBDGPU : public TwoStepBD
    {
    public:
        //! Constructs the integration method and associates it with the system
        TwoStepBDGPU(std::shared_ptr<SystemDefinition> sysdef,
                     std::shared_ptr<ParticleGroup> group,
                     std::shared_ptr<Variant> T,
                     unsigned int seed,
                     bool use_lambda,
                     Scalar lambda,
                     bool noiseless_t,
                     bool noiseless_r);
        ~TwoStepBDGPU();

        //! Performs the first step of the integration
        virtual void integrateStepOne(unsigned int timestep);

        //! Performs the second step of the integration
        virtual void integrateStepTwo(unsigned int timestep);
        
        //void readData();
        std::fstream pitsFile;
        std::string name = "/home/ilya/MD/TestHD/data.txt";
        
        int nMax1 = 0;                                   // realinoe 4islo zone
        const int SIZE = 1000;                           // max 4islo zone
        double *x_pitCPU;                                // massiv dl9 cpu
        double *y_pitCPU;
        double *x_pit;                                   // massivi dl9 gpu
        double *y_pit;
        double mySigma1, epsilon1;                       // zone coefs
        
    protected:
        unsigned int m_block_size;               //!< block size
    };

//! Exports the TwoStepBDGPU class to python
void export_TwoStepBDGPU(pybind11::module& m);

#endif // #ifndef __TWO_STEP_BD_GPU_H__
