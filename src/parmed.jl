using PyCall

function MMSystem(parm7::AbstractString, rst7::AbstractString; cutoff=nothing)
    amb = pyimport("parmed.amber")
    parm = amb.AmberParm(parm7, rst7)

    # hard code Ewald parameters for now
    alpha = 0.394670228244391 * 10
    kmax = 20
    spline_order = 6

    if isnothing(parm.box)
        box = [0., 0., 0.]
        recip = nothing
    elseif parm.box[4:6] != [90., 90., 90.]
        error("Only orthogonal box is supported")
    else
        box = parm.box[1:3] / 10.0
        gridpoints = Tuple(floor.(Int, parm.box[1:3]))
        # ewald = EwaldRecip(alpha, length(parm.atoms), kmax, box)
        recip = PMERecip(alpha, spline_order, gridpoints, JuMD.KE, Threads.nthreads())
    end

    positions = [(x._value ./ 10.0) for x in parm.positions]
    masses = [atom.mass for atom in parm.atoms]
    atomic_numbers = [atom.atomic_number for atom in parm.atoms]
    nonbonded_exclusion = [Set{Int}() for _ in parm.atoms]

    bonds = HarmonicBondForce()
    for bond in parm.bonds
        index1 = bond.atom1.idx + 1
        index2 = bond.atom2.idx + 1
        push!(bonds.indices, (index1, index2))
        push!(bonds.length, bond.type.req / 10.0)
        push!(bonds.k, bond.type.k * 2 * 4.184 * 100)
        push!(nonbonded_exclusion[index1], index2)
        push!(nonbonded_exclusion[index2], index1)
    end

    angles = HarmonicAngleForce()
    for angle in parm.angles
        index1 = angle.atom1.idx + 1
        index2 = angle.atom2.idx + 1
        index3 = angle.atom3.idx + 1
        push!(angles.indices, (index1, index2, index3))
        push!(angles.angle, angle.type.theteq / 180.0 * pi)
        push!(angles.k, angle.type.k * 4.184 * 2)
        push!(nonbonded_exclusion[index1], index3)
        push!(nonbonded_exclusion[index3], index1)
    end

    dihedrals = PeriodicTorsionForce()
    for dihedral in parm.dihedrals
        index1 = dihedral.atom1.idx + 1
        index2 = dihedral.atom2.idx + 1
        index3 = dihedral.atom3.idx + 1
        index4 = dihedral.atom4.idx + 1
        push!(dihedrals.indices, (index1, index2, index3, index4))
        push!(dihedrals.periodicity, dihedral.type.per)
        push!(dihedrals.phase, dihedral.type.phase / 180.0 * pi)
        push!(dihedrals.k, dihedral.type.phi_k * 4.184)
        if !dihedral.improper
            push!(nonbonded_exclusion[index1], index4)
            push!(nonbonded_exclusion[index4], index1)
        end
    end

    vdw14 = LennardJonesExceptionForce()
    elec14 = CoulombExceptionForce()
    indice_set = Set()
    for dihedral in parm.dihedrals
        if !dihedral.improper
            index1 = dihedral.atom1.idx + 1
            index2 = dihedral.atom4.idx + 1
            if (index1, index2) ∉ indice_set
                push!(indice_set, (index1, index2))

                sigma1 = parm.atoms[index1].sigma
                sigma2 = parm.atoms[index2].sigma
                epsilon1 = parm.atoms[index1].epsilon
                epsilon2 = parm.atoms[index2].epsilon
                sigma = (sigma1 + sigma2) / 2
                epsilon = 4 * sqrt(epsilon1 * epsilon2) / dihedral.type.scnb
                push!(vdw14.indices, (index1, index2))
                push!(vdw14.sigma, sigma / 10.0)
                push!(vdw14.epsilon, epsilon * 4.184)

                q1 = parm.atoms[index1].charge
                q2 = parm.atoms[index2].charge
                charge_prod = q1 * q2 / dihedral.type.scee
                push!(elec14.indices, (index1, index2))
                push!(elec14.charge_prod, charge_prod)
            end
        end
    end

    nonbonded_exclusion = collect.(nonbonded_exclusion)

    vdw = LennardJonesForce(exclusion=nonbonded_exclusion, cutoff=cutoff)
    elec = CoulombForce(exclusion=nonbonded_exclusion, cutoff=cutoff, recip=recip)
    for atom in parm.atoms
        push!(vdw.sigma, atom.sigma / 10.0 / 2)
        push!(vdw.epsilon, 2 * sqrt(atom.epsilon * 4.184))
        push!(elec.charges, atom.charge)
    end

    force_groups = ForceGroups(groups=(bonds=bonds, angles=angles, dihedrals=dihedrals, vdw=vdw, vdw14=vdw14, elec=elec, elec14=elec14))

    MMSystem(box, positions, masses, atomic_numbers, force_groups)
end
