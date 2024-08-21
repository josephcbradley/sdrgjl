const σ_x = Hermitian([0.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 0.0 + 0.0im]) |> SparseMatrixCSC

const σ_y = Hermitian([0.0 + 0.0im 0.0 + -1.0im; 0.0 + 1.0im 0.0 + 0.0im]) |> SparseMatrixCSC

const σ_z = Hermitian([1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im -1.0 + 0.0im]) |> SparseMatrixCSC

## Assume h_bar = 1.0

const Sx = 0.5 * σ_x

const Sy = 0.5 * σ_y

const Sz = 0.5 * σ_z

const up = @SVector [1.0, 0.0]

const down = @SVector [0.0, 1.0]

const upup = kron(up, up)

const updown = kron(up, down)

const downup = kron(down, up)

const downdown = kron(down, down)
