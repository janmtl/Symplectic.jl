module Symplectic

  import Base.size
  import Base.getindex
  import Base.convert
  import Base.ndims
  import Base.size
  import Base.repeat

	import Base.*
	import Base./
	import Base.==

  # New symplectic matrix type
  type Symp
      x::Array{Bool, 2}
      z::Array{Bool, 2}
  end

  # Overload the basic matrix methods
  function size(s::Symp)
    return size(s.x)
  end

  function getindex(s::Symp, i, j)
    return Symp(convert(Array{Bool, 2}, s.x[i,j]), convert(Array{Bool, 2}, s.z[i,j]))
  end

  function convert(::Type{Array{Bool, 2}}, s::Bool)
    return fill(s, (1, 1))::Array{Bool, 2}
  end

  function ndims(s::Symp)
    return 2
  end

  function size(s::Symp, t::Int64)
    return size(s.x::Array{Bool, 2}, t::Int64)
  end

  function repeat(s::Symp; inner=[1 1], outer=[1 1])
    return Symp(repeat(s.x, inner=inner, outer=outer), repeat(s.z, inner=inner, outer=outer))
  end

  # Commutation
  function *(a::Symp, b::Symp)
    n_a = size(b)[1];
    n_b = size(a)[1];
    r_a = repeat(a, inner=[n_a, 1]);
    r_b = repeat(b, outer=[n_b, 1]);
    res = reducedim($, r_a.x & r_b.z, [2], false) $ reducedim($, r_a.z & r_b.x, [2], false);
    return reshape(res, (n_a, n_b))
  end

  # Anti-commutation
  function /(a::Symp, b::Symp)
    return ~(*(a::Symp, b::Symp))
  end

  # Equivalence
  function ==(a::Symp, b::Symp)
    eq_x = reducedim(&, a.x .== b.x, [2], true)
    eq_z = reducedim(&, a.z .== b.z, [2], true)
    return (eq_x & eq_z)[1]
  end

  export Symp, size, getindex, convert, ndims, size, repeat, *, /
end # module
