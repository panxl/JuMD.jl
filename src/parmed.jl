using PyCall

function MMSystem(parm7::AbstractString, rst7::AbstractString)
    amb = pyimport("parmed.amber")
    parm = amb.AmberParm(parm7, rst7)

    positions = Vector{SVector{3, Float64}}(parm.positions._value)
    velocities = similar(positions)
    fill!(velocities, zeros(eltype(velocities)))
    forces = similar(positions)
    fill!(forces, zeros(eltype(forces)))
    masses = [atom.mass for atom in parm.atoms]
    atomic_numbers = [atom.atomic_number for atom in parm.atoms]
    nonbonded_exslusion = [Set{Int}() for _ in parm.atoms]

    force_group_bonds = ForceGroup(HarmonicBondForce[])
    for bond in parm.bonds
        push!(force_group_bonds.forces, HarmonicBondForce((bond.atom1.idx+1, bond.atom2.idx+1), bond.type.req, bond.type.k*2))
        push!(nonbonded_exslusion[bond.atom1.idx+1], bond.atom2.idx+1)
        push!(nonbonded_exslusion[bond.atom2.idx+1], bond.atom1.idx+1)
    end

    force_group_angles = ForceGroup(HarmonicAngleForce[])
    for angle in parm.angles
        push!(force_group_angles.forces, HarmonicAngleForce((angle.atom1.idx+1, angle.atom2.idx+1, angle.atom3.idx+1), angle.type.theteq/180.0*pi, angle.type.k*2))
        push!(nonbonded_exslusion[angle.atom1.idx+1], angle.atom3.idx+1)
        push!(nonbonded_exslusion[angle.atom3.idx+1], angle.atom1.idx+1)
    end

    force_group_dihedrals = ForceGroup(PeriodicTorsionForce[])
    for dihedral in parm.dihedrals
        push!(force_group_dihedrals.forces, PeriodicTorsionForce((dihedral.atom1.idx+1, dihedral.atom2.idx+1, dihedral.atom3.idx+1, dihedral.atom4.idx+1), dihedral.type.per, dihedral.type.phase/180.0*pi, dihedral.type.phi_k))
        if !dihedral.improper
            push!(nonbonded_exslusion[dihedral.atom1.idx+1], dihedral.atom4.idx+1)
            push!(nonbonded_exslusion[dihedral.atom4.idx+1], dihedral.atom1.idx+1)
        end
    end

    vdw_force = LennardJonesForce([atom.sigma/2 for atom in parm.atoms], [2*sqrt(atom.epsilon) for atom in parm.atoms], nonbonded_exslusion)
    force_group_vdw = ForceGroup(LennardJonesForce[vdw_force])    

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
            push!(vdw14_force, LennardJonesExceptionForce((index1, index2), sigma, epsilon))
        end
    end
    force_group_vdw14 = ForceGroup(collect(vdw14_force))

    elec_force = CoulombForce([atom.charge for atom in parm.atoms], nonbonded_exslusion)
    force_group_elec = ForceGroup(CoulombForce[elec_force])    

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

    MMSystem(positions, velocities, forces, masses, atomic_numbers, force_groups)
end
