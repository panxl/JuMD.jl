function shift(r, rcut)
    r_rcut = r / rcut
    r_rcut² = r_rcut^2
    s = (1 - r_rcut²)^2
    dsdr = -4 * (1 - r_rcut²) * r_rcut / rcut
    return s, dsdr
end
