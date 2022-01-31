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

    force_group_bonds = ForceGroup(HarmonicBondForce[])
    for bond in parm.bonds
        push!(force_group_bonds.forces, HarmonicBondForce((bond.atom1.idx+1, bond.atom2.idx+1), bond.type.req, bond.type.k*2))
    end

    force_group_angles = ForceGroup(HarmonicAngleForce[])
    for angle in parm.angles
        push!(force_group_angles.forces, HarmonicAngleForce((angle.atom1.idx+1, angle.atom2.idx+1, angle.atom3.idx+1), angle.type.theteq/180.0*pi, angle.type.k*2))
    end

    force_group_dihedrals = ForceGroup(PeriodicTorsionForce[])
    for dihedral in parm.dihedrals
        push!(force_group_dihedrals.forces, PeriodicTorsionForce((dihedral.atom1.idx+1, dihedral.atom2.idx+1, dihedral.atom3.idx+1, dihedral.atom4.idx+1), dihedral.type.per, dihedral.type.phase/180.0*pi, dihedral.type.phi_k))
    end

    force_group = ForceGroups((bonds=force_group_bonds, angles=force_group_angles, dihedrals=force_group_dihedrals))

    MMSystem(positions, velocities, forces, masses, atomic_numbers, force_group)
end