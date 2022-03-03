module helpme

using helpme_jll

export PMEInstance, helpme_createD, helpme_destroyD
export helpme_setupD, helpme_setup_compressedD, helpme_set_lattice_vectorsD
export helpme_compute_E_recD, helpme_compute_EF_recD
export helpme_compute_EFV_recD, helpme_compute_P_recD
export XAligned, ShapeMatrix

@enum LatticeType XAligned=1 ShapeMatrix=2

struct PMEInstance end

function helpme_createD()
    @ccall libhelpme.helpme_createD()::Ptr{PMEInstance}
end

function helpme_destroyD(pme::Ptr{PMEInstance})
    @ccall libhelpme.helpme_destroyD(pme::Ptr{PMEInstance})::Cvoid
end

function helpme_setupD(pme::Ptr{PMEInstance}, rPower, kappa, splineOrder, aDim, bDim, cDim, scaleFactor, nThreads)
    @ccall libhelpme.helpme_setupD(
        pme::Ptr{PMEInstance},
        rPower::Cshort,
        kappa::Cdouble,
        splineOrder::Cint,
        aDim::Cint,
        bDim::Cint,
        cDim::Cint,
        scaleFactor::Cdouble,
        nThreads::Cint,
    )::Cvoid
end

function helpme_setup_compressedD(pme::Ptr{PMEInstance}, rPower, kappa, splineOrder, aDim, bDim, cDim, maxKA, maxKB, maxKC, scaleFactor, nThreads)
    @ccall libhelpme.helpme_setup_compressedD(
        pme::Ptr{PMEInstance},
        rPower::Cshort,
        kappa::Cdouble,
        splineOrder::Cint,
        aDim::Cint,
        bDim::Cint,
        cDim::Cint,
        maxKA::Cint,
        maxKB::Cint,
        maxKC::Cint,
        scaleFactor::Cdouble,
        nThreads::Cint,
    )::Cvoid
end

function helpme_set_lattice_vectorsD(pme::Ptr{PMEInstance}, A, B, C, alpha, beta, gamma, LatticeType)
    @ccall libhelpme.helpme_set_lattice_vectorsD(
        pme::Ptr{PMEInstance},
        A::Cdouble,
        B::Cdouble,
        C::Cdouble,
        alpha::Cdouble,
        beta::Cdouble,
        gamma::Cdouble,
        LatticeType::Cint,
    )::Cvoid
end

function helpme_compute_E_recD(pme::Ptr{PMEInstance}, nAtoms, parameterAngMom, parameters, coordinates)
    @ccall libhelpme.helpme_compute_E_recD(
        pme::Ptr{PMEInstance},
        nAtoms::Csize_t,
        parameterAngMom::Cint,
        parameters::Ref{Cdouble},
        coordinates::Ref{Cdouble},
    )::Cdouble
end

function helpme_compute_EF_recD(pme::Ptr{PMEInstance}, nAtoms, parameterAngMom, parameters, coordinates, forces)
    @ccall libhelpme.helpme_compute_EF_recD(
        pme::Ptr{PMEInstance},
        nAtoms::Csize_t,
        parameterAngMom::Cint,
        parameters::Ref{Cdouble},
        coordinates::Ref{Cdouble},
        forces::Ref{Cdouble},
    )::Cdouble
end

function helpme_compute_EFV_recD(pme::Ptr{PMEInstance}, nAtoms, parameterAngMom, parameters, coordinates, forces, virial)
    @ccall libhelpme.helpme_compute_EFV_recD(
        pme::Ptr{PMEInstance},
        nAtoms::Csize_t,
        parameterAngMom::Cint,
        parameters::Ref{Cdouble},
        coordinates::Ref{Cdouble},
        forces::Ref{Cdouble},
        virial::Ref{Cdouble},
    )::Cdouble
end

function helpme_compute_P_recD(pme::Ptr{PMEInstance}, nAtoms, parameterAngMom, parameters, coordinates, nGridPoints, gridPoints, derivativeLevel, potential)
    @ccall libhelpme.helpme_compute_P_recD(
        pme::Ptr{PMEInstance},
        nAtoms::Csize_t,
        parameterAngMom::Cint,
        parameters::Ref{Cdouble},
        coordinates::Ref{Cdouble},
        nGridPoints::Csize_t,
        gridPoints::Ref{Cdouble},
        derivativeLevel::Cint,
        potential::Ref{Cdouble},
    )::Cvoid
end

end # module
