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

    positions = [SVector(x._value ./ 10.0) for x in parm.positions]
    masses = [atom.mass for atom in parm.atoms]
    atomic_numbers = [atom.atomic_number for atom in parm.atoms]
    nonbonded_exslusion = [Set{Int}() for _ in parm.atoms]

    force_group_bonds = ForceGroup(HarmonicBondForce[])
    for bond in parm.bonds
        push!(force_group_bonds.forces, HarmonicBondForce((bond.atom1.idx+1, bond.atom2.idx+1), bond.type.req/10.0, bond.type.k*2*4.184*100))
        push!(nonbonded_exslusion[bond.atom1.idx+1], bond.atom2.idx+1)
        push!(nonbonded_exslusion[bond.atom2.idx+1], bond.atom1.idx+1)
    end

    force_group_angles = ForceGroup(HarmonicAngleForce[])
    for angle in parm.angles
        push!(force_group_angles.forces, HarmonicAngleForce((angle.atom1.idx+1, angle.atom2.idx+1, angle.atom3.idx+1), angle.type.theteq/180.0*pi, angle.type.k*4.184*2))
        push!(nonbonded_exslusion[angle.atom1.idx+1], angle.atom3.idx+1)
        push!(nonbonded_exslusion[angle.atom3.idx+1], angle.atom1.idx+1)
    end

    force_group_dihedrals = ForceGroup(PeriodicTorsionForce[])
    for dihedral in parm.dihedrals
        push!(force_group_dihedrals.forces, PeriodicTorsionForce((dihedral.atom1.idx+1, dihedral.atom2.idx+1, dihedral.atom3.idx+1, dihedral.atom4.idx+1), dihedral.type.per, dihedral.type.phase/180.0*pi, dihedral.type.phi_k*4.184))
        if !dihedral.improper
            push!(nonbonded_exslusion[dihedral.atom1.idx+1], dihedral.atom4.idx+1)
            push!(nonbonded_exslusion[dihedral.atom4.idx+1], dihedral.atom1.idx+1)
        end
    end

    nonbonded_exslusion = collect.(nonbonded_exslusion)

    vdw_force = LennardJonesForce([atom.sigma/10.0/2 for atom in parm.atoms], [2*sqrt(atom.epsilon*4.184) for atom in parm.atoms], nonbonded_exslusion)
    force_group_vdw = ForceGroup(vdw_force)

    vdw14_force = Set{LennardJonesExceptionForce}()
    for dihedral in parm.dihedrals
        if !dihedral.improper
            index1 = dihedral.atom1.idx + 1
            index2 = dihedral.atom4.idx + 1
            sigma1 = parm.atoms[index1].sigma
            sigma2 = parm.atoms[index2].sigma
            epsilon1 = parm.atoms[index1].epsilon
            epsilon2 = parm.atoms[index2].epsilon
            sigma = (sigma1 + sigma2) / 2
            epsilon = 4 * sqrt(epsilon1 * epsilon2) / dihedral.type.scnb
            push!(vdw14_force, LennardJonesExceptionForce((index1, index2), sigma/10.0, epsilon*4.184))
        end
    end
    force_group_vdw14 = ForceGroup(collect(vdw14_force))

    elec_force = CoulombForce([atom.charge for atom in parm.atoms], nonbonded_exslusion, recip)
    force_group_elec = ForceGroup(elec_force)

    elec14_force = Set{CoulombExceptionForce}()
    for dihedral in parm.dihedrals
        if !dihedral.improper
            index1 = dihedral.atom1.idx + 1
            index2 = dihedral.atom4.idx + 1
            q1 = parm.atoms[index1].charge
            q2 = parm.atoms[index2].charge
            push!(elec14_force, CoulombExceptionForce((index1, index2), (q1, q2), 1 / dihedral.type.scee))
        end
    end
    force_group_elec14 = ForceGroup(collect(elec14_force))

    force_groups = ForceGroups((bonds=force_group_bonds, angles=force_group_angles, dihedrals=force_group_dihedrals, vdw=force_group_vdw, vdw14=force_group_vdw14, elec=force_group_elec, elec14=force_group_elec14))

    MMSystem(box, positions, masses, atomic_numbers, force_groups, cutoff=cutoff)
end
