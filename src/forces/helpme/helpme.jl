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
    ccall((:helpme_createD, libhelpme), Ptr{PMEInstance}, ())
end

function helpme_destroyD(pme::Ptr{PMEInstance})
    ccall((:helpme_destroyD, libhelpme), Cvoid, (Ptr{PMEInstance},), pme,)
end

function helpme_setupD(pme::Ptr{PMEInstance}, rPower, kappa, splineOrder, aDim, bDim, cDim, scaleFactor, nThreads)
    ccall(
        (:helpme_setupD, libhelpme), Cvoid,
        (Ptr{PMEInstance}, Cshort, Cdouble, Cint, Cint, Cint, Cint, Cdouble, Cint),
        pme, rPower, kappa, splineOrder, aDim, bDim, cDim, scaleFactor, nThreads,
        )
end

function helpme_setup_compressedD(pme::Ptr{PMEInstance}, rPower, kappa, splineOrder, aDim, bDim, cDim, maxKA, maxKB, maxKC, scaleFactor, nThreads)
    ccall(
        (:helpme_setup_compressedD, libhelpme), Cvoid,
        (Ptr{PMEInstance}, Cshort, Cdouble, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cdouble, Cint),
        pme, rPower, kappa, splineOrder, aDim, bDim, cDim, maxKA, maxKB, maxKC, scaleFactor, nThreads,
        )
end

function helpme_set_lattice_vectorsD(pme::Ptr{PMEInstance}, A, B, C, alpha, beta, gamma, LatticeType)
    ccall(
        (:helpme_set_lattice_vectorsD, libhelpme), Cvoid,
        (Ptr{PMEInstance}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cint),
        pme, A, B, C, alpha, beta, gamma, LatticeType,
        )
end

function helpme_compute_E_recD(pme::Ptr{PMEInstance}, nAtoms, parameterAngMom, parameters, coordinates)
    ccall(
        (:helpme_compute_E_recD, libhelpme), Cdouble,
        (Ptr{PMEInstance}, Csize_t, Cint, Ref{Cdouble}, Ref{Cdouble}),
        pme, nAtoms, parameterAngMom, parameters, coordinates,
        )
end

function helpme_compute_EF_recD(pme::Ptr{PMEInstance}, nAtoms, parameterAngMom, parameters, coordinates, forces)
    ccall(
        (:helpme_compute_EF_recD, libhelpme), Cdouble,
        (Ptr{PMEInstance}, Csize_t, Cint, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}),
        pme, nAtoms, parameterAngMom, parameters, coordinates, forces,
        )
end

function helpme_compute_EFV_recD(pme::Ptr{PMEInstance}, nAtoms, parameterAngMom, parameters, coordinates, forces, virial)
    ccall(
        (:helpme_compute_EFV_recD, libhelpme), Cdouble,
        (Ptr{PMEInstance}, Csize_t, Cint, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}),
        pme, nAtoms, parameterAngMom, parameters, coordinates, forces, virial,
        )
end

function helpme_compute_P_recD(pme::Ptr{PMEInstance}, nAtoms, parameterAngMom, parameters, coordinates, nGridPoints, gridPoints, derivativeLevel, potential)
    ccall(
        (:helpme_compute_P_recD, libhelpme), Cvoid,
        (Ptr{PMEInstance}, Csize_t, Cint, Ref{Cdouble}, Ref{Cdouble}, Csize_t, Ref{Cdouble}, Cint, Ref{Cdouble}),
        pme, nAtoms, parameterAngMom, parameters, coordinates, nGridPoints, gridPoints, derivativeLevel, potential,
        )
end

end # module
