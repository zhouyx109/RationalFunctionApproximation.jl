# Certify poles and uniform errors of the stored interval fits.
# Treat each stored Float64 coefficient as its exact binary rational.
include("run.jl")

const QQ=Rational{BigInt}
const EQ=Complex{QQ}

exact(z::Number)=EQ(QQ(real(z)),QQ(imag(z)))
degree(p)=length(p)-1
function trim(p)
    p=copy(p)
    while length(p)>1 && iszero(last(p));pop!(p);end
    p
end
iszero_poly(p)=all(iszero,p)
psub(a,b)=padd(a,-b)
function padd(a,b)
    c=fill(zero(promote_type(eltype(a),eltype(b))),max(length(a),length(b)))
    c[1:length(a)].+=a;c[1:length(b)].+=b;trim(c)
end
function pmul(a,b)
    c=fill(zero(promote_type(eltype(a),eltype(b))),length(a)+length(b)-1)
    for i in eachindex(a),j in eachindex(b);c[i+j-1]+=a[i]*b[j];end
    trim(c)
end
function pdivrem(a,b)
    iszero_poly(b) && error("zero polynomial divisor")
    r=trim(a);q=fill(zero(eltype(r)),max(length(r)-length(b)+1,1))
    while !iszero_poly(r) && degree(r)>=degree(b)
        k=degree(r)-degree(b);v=last(r)/last(b);q[k+1]+=v
        for j in eachindex(b);r[k+j]-=v*b[j];end
        r=trim(r)
    end
    trim(q),r
end
function pgcd(a,b)
    a,b=trim(a),trim(b)
    while !iszero_poly(b);a,b=b,last(pdivrem(a,b));end
    iszero_poly(a) ? a : trim(a./last(a))
end
function pexactdiv(a,b)
    q,r=pdivrem(a,b);iszero_poly(r)||error("nonexact polynomial division");q
end
peval(p,z)=foldr((a,b)->a+z*b,p)
pderiv(p)=length(p)==1 ? [zero(first(p))] : [j*p[j+1] for j in 1:length(p)-1]

"Exact stored-coefficient model matching finite-point evaluation semantics."
function exact_blocks(r;expanded=false)
    d=r isa RFA.QuadraticThiele ? getfield(r,:data) : r
    if d isa Thiele
        b=[[exact(w)] for w in weights(d)]
        a=[[exact(-z),one(EQ)] for z in nodes(d)]
    else
        b=[exact.(v) for v in d.bpolys]
        a=Vector{EQ}[]
        for (v,pair) in zip(d.apolys,d.pairs)
            if expanded
                push!(a,exact.(v))
            elseif length(v)==2
                push!(a,exact(v[2]).*EQ[-exact(pair[1]),one(EQ)])
            else
                push!(a,exact(v[3]).*pmul(EQ[-exact(pair[1]),one(EQ)],EQ[-exact(pair[2]),one(EQ)]))
            end
        end
    end
    b,a
end

function chart_continuants(r, coefficient, sum_poly, multiply_poly)
    groups=r.chart.groups
    one_coefficient=coefficient(1.0)
    linear(x)=[coefficient(-x),one_coefficient]
    complement(g)=foldr((x,p) -> sum_poly([coefficient(x)],multiply_poly(linear(g.center),p)),
        g.complement;init=[zero(one_coefficient)])
    g=last(groups)
    p,q=isempty(g.gamma) ? (coefficient.(g.head),[one_coefficient]) :
        (complement(g),coefficient.(g.gamma))
    for k in length(groups)-1:-1:1
        g,pair=groups[k],r.pairs[2k]
        a=linear(pair[1])
        pair[1]==pair[2] || (a=multiply_poly(a,linear(pair[2])))
        rp=multiply_poly([coefficient(g.ratio)],p)
        aq=multiply_poly(a,q)
        p,q=sum_poly(multiply_poly(complement(g),rp),multiply_poly(coefficient.(g.head),aq)),
            sum_poly(multiply_poly(coefficient.(g.gamma),rp),aq)
    end
    trim(p),trim(q)
end

function exact_chart_determinants(r)
    factors=Vector{EQ}[]
    for k in 1:length(r.chart.groups)-1
        g,pair=r.chart.groups[k],r.pairs[2k]
        push!(factors,EQ[-exact(pair[1]),1])
        pair[1]==pair[2] || push!(factors,EQ[-exact(pair[2]),1])
        c=foldr((x,p) -> padd(EQ[exact(x)],pmul(EQ[-exact(g.center),1],p)),
            g.complement;init=EQ[0])
        push!(factors,psub(c,pmul(exact.(g.head),exact.(g.gamma))))
    end
    factors
end

function exact_continuants(r;expanded=false)
    if r isa RFA.QuadraticThiele && r.chart !== nothing
        p,q=chart_continuants(r,exact,padd,pmul)
        return p,q,exact_chart_determinants(r)
    end
    b,a=exact_blocks(r;expanded)
    isempty(b) && return EQ[0],EQ[1],Vector{EQ}[]
    po,qo=EQ[1],EQ[0];p,q=copy(first(b)),EQ[1]
    for j in 2:length(b)
        p,q,po,qo=padd(pmul(b[j],p),pmul(a[j-1],po)),
            padd(pmul(b[j],q),pmul(a[j-1],qo)),p,q
    end
    trim(p),trim(q),a[1:max(length(b)-1,0)]
end

"GCD factors are sought only in the exact stored continuant determinant product."
function reduce_exact(p,q,a)
    g=EQ[1];factors=Vector{EQ}[]
    for aj in a
        if iszero_poly(aj)
            h=pgcd(p,q)
            return pexactdiv(p,h),pexactdiv(q,h),pmul(g,h),vcat(factors,[h])
        end
        h=pgcd(aj,pgcd(last(pdivrem(p,aj)),last(pdivrem(q,aj))))
        if degree(h)>0
            p=pexactdiv(p,h);q=pexactdiv(q,h);g=pmul(g,h);push!(factors,h)
        end
    end
    p,q,g,factors
end
function lift_poly(p)
    a=fill(zero(eltype(p)),2length(p)-1)
    for j in eachindex(p);a[2j-1]=p[j];end
    trim(a)
end
function exact_model(r,kind;expanded=false)
    pu,qu,a=exact_continuants(r;expanded)
    p,q,g,factors=reduce_exact(pu,qu,a)
    if kind in (:even,:odd)
        p,q=lift_poly(p),lift_poly(q)
        g=lift_poly(g)
        kind===:odd && (p=pmul(EQ[0,1],p))
        while length(p)>1 && length(q)>1 && iszero(first(p)) && iszero(first(q))
            p=p[2:end];q=q[2:end];g=pmul(g,EQ[0,1])
        end
        pu,qu=pmul(p,g),pmul(q,g)
    end
    (;p=trim(p),q=trim(q),g=trim(g),pu,qu,factors)
end

# Rotate the fitting axis onto the real line before root counting.
function axis_polynomials(m,phase)
    p=[m.p[j]*phase^(j-1) for j in eachindex(m.p)]
    q=[m.q[j]*phase^(j-1) for j in eachindex(m.q)]
    nonzero_constant=!iszero(first(q))
    scale=nonzero_constant ? first(q) : last(q)
    p,q=p./scale,q./scale
    cert=nonzero_constant && all(j->iseven(j-1) ? isreal(q[j])&&real(q[j])>=0 : iszero(q[j]),eachindex(q))
    p,q,cert
end

# Integer-scaled Bernstein coefficients on [0,1]. Sign variation bounds
# the number of roots in the open interval by Descartes after u=t/(1+t).
function bernstein_integers(c)
    n=length(c)-1
    a=[c[j+1]/binomial(big(n),j) for j in 0:n]
    D=foldl(lcm,denominator.(a);init=BigInt(1))
    A=BigInt[numerator(v)*(D÷denominator(v)) for v in a]
    BigInt[sum(binomial(big(i),j)*A[j+1] for j in 0:i) for i in 0:n]
end
function variations(v)
    old=0;n=0
    for x in v
        s=sign(x);iszero(s)&&continue
        !iszero(old)&&s!=old&&(n+=1);old=s
    end
    n
end
function split_bernstein(B)
    n=length(B)-1;t=copy(B);left=copy(B);right=copy(B)
    left[1]=t[1]<<n;right[n+1]=t[n+1]<<n
    for j in 1:n
        t=[t[k]+t[k+1] for k in 1:length(t)-1]
        left[j+1]=t[1]<<(n-j);right[n-j+1]=t[end]<<(n-j)
    end
    left,right
end
function exact_axis_count(q;maxdepth=512,maxnodes=20000)
    even=all(j->iseven(j-1)||iszero(q[j]),eachindex(q))
    if !even || !all(isreal,q)
        return (;status="unsupported_axis_polynomial",rootcount=-1,nodes=0,maxdepth=0,initial_variations=-1,
                brackets=Tuple{QQ,QQ}[],endpoints=QQ[])
    end
    count_unit(real.(q[1:2:end]);maxdepth,maxnodes)
end

function count_unit(c;maxdepth=512,maxnodes=20000,maxseconds=15.0)
    B=bernstein_integers(c);started=time()
    initial=variations(B);stack=[(QQ(0),QQ(1),B,0)];seen=0;depth=0
    endpoints=Set{QQ}();brackets=Tuple{QQ,QQ}[]
    while !isempty(stack)
        a,b,v,l=pop!(stack);seen+=1;depth=max(depth,l)
        iszero(first(v))&&push!(endpoints,a);iszero(last(v))&&push!(endpoints,b)
        count=variations(v)
        count==0&&continue
        if count==1
            push!(brackets,(a,b));continue
        end
        if l>=maxdepth || seen>=maxnodes || time()-started>=maxseconds
            return (;status="unresolved_variations",rootcount=-1,nodes=seen,maxdepth=depth,initial_variations=initial,
                    brackets,endpoints=collect(endpoints))
        end
        L,R=split_bernstein(v);mid=(a+b)/2
        push!(stack,(a,mid,L,l+1));push!(stack,(mid,b,R,l+1))
    end
    (;status="exact_bernstein_descartes",rootcount=length(brackets)+length(endpoints),nodes=seen,maxdepth=depth,
        initial_variations=initial,brackets,endpoints=collect(endpoints))
end

const CHECKPRIMES=(65537,65539,65543)
function modular_coprime(p,q)
    for ell in CHECKPRIMES
        # Dividing by the constant coefficient can introduce non-dyadic denominators.
        any(c->iszero(mod(denominator(real(c)),ell)),p) && continue
        any(c->iszero(mod(denominator(real(c)),ell)),q) && continue
        modc(c)=Int[mod(numerator(real(x))*invmod(mod(denominator(real(x)),ell),ell),ell) for x in c]
        a,b=modc(p),modc(q)
        (iszero(last(a))||iszero(last(b)))&&continue
        while !(length(b)==1 && iszero(first(b)))
            v=invmod(last(b),ell)
            while length(a)>=length(b) && !(length(a)==1 && iszero(first(a)))
                d=length(a)-length(b);t=mod(last(a)*v,ell)
                for j in eachindex(b);a[j+d]=mod(a[j+d]-t*b[j],ell);end
                while length(a)>1&&iszero(last(a));pop!(a);end
            end
            a,b=b,a
        end
        length(a)==1&&!iszero(first(a))&&return ell
    end
    0
end

function strict_positive(c;maxseconds=15.0,maxnodes=20000,maxdepth=512)
    v=bernstein_integers(c);stack=[(QQ(0),QQ(1),v,0)];started=time();n=0;dep=0
    while !isempty(stack)
        a,b,w,l=pop!(stack);n+=1;dep=max(dep,l)
        if first(w)<=0 || last(w)<=0
            x=first(w)<=0 ? a : b
            @assert peval(c,x)<=0
            return (;status="nonpositive_witness",n,dep,witness=string(x),seconds=time()-started)
        end
        variations(w)==0&&continue
        if n>=maxnodes||l>=maxdepth||time()-started>maxseconds
            return (;status="unresolved_budget",n,dep,witness="",seconds=time()-started)
        end
        L,R=split_bernstein(w);m=(a+b)/2
        push!(stack,(a,m,L,l+1));push!(stack,(m,b,R,l+1))
    end
    (;status="strictly_positive",n,dep,witness="",seconds=time()-started)
end

const CERT_THRESHOLDS=(QQ(1,10^11),QQ(3,2*10^11),QQ(2,10^11),QQ(4,10^11),QQ(5,10^11),QQ(1,10^10))
const SINE_CERT_THRESHOLDS=(QQ(1,10^13),QQ(2,10^13),QQ(5,10^13),QQ(1,10^12),CERT_THRESHOLDS...)
# The Taylor polynomial through degree 18 has degree 17. On [-1,1],
# the remainder has absolute value at most 1/19!.
const SINE_TAYLOR=QQ[iseven(j) ? 0 : QQ((-1)^((j-1)÷2),factorial(big(j))) for j in 0:17]
const SINE_TAYLOR_REMAINDER=QQ(1,factorial(big(19)))

function save_cert(outdir,name,rows)
    isempty(rows)&&return
    path=joinpath(outdir,name)
    open(path*".tmp","w") do io
        println(io,join(keys(first(rows)),","))
        for r in rows
            vals=map(values(r)) do val
                str=string(val)
                occursin(r"[,\"\n]",str) ? "\""*replace(str,"\""=>"\"\"")*"\"" : str
            end
            println(io,join(vals,","))
        end
    end
    mv(path*".tmp",path;force=true)
end

# Use a backward recurrence to check the exact forward continuants.
# Lifted models evaluate the scalar fraction at z² first.
function exact_chart_pair(r,z)
    complement(g)=evalpoly(z-exact(g.center),exact.(g.complement))
    groups=r.chart.groups
    g=last(groups)
    n,d=isempty(g.gamma) ? (peval(exact.(g.head),z),one(EQ)) :
        (complement(g),peval(exact.(g.gamma),z))
    for k in length(groups)-1:-1:1
        g,pair=groups[k],r.pairs[2k]
        a=z-exact(pair[1])
        pair[1]==pair[2] || (a*=z-exact(pair[2]))
        u,v=exact(g.ratio)*n,a*d
        n,d=complement(g)*u+peval(exact.(g.head),z)*v,
            peval(exact.(g.gamma),z)*u+v
    end
    n,d
end

function check_model(r,kind,p,q,phase)
    chart=r isa RFA.QuadraticThiele && r.chart !== nothing
    bs,as=chart ? (nothing,nothing) : exact_blocks(r)
    for x in (QQ(0),QQ(1,17),QQ(-2,5),QQ(1))
        z=phase*x
        w=kind in (:even,:odd) ? z*z : z
        if chart
            n,d=exact_chart_pair(r,w)
        else
            n,d=peval(last(bs),w),one(EQ)
            for j in length(bs)-1:-1:1
                n,d=peval(bs[j],w)*n+peval(as[j],w)*d,n
            end
        end
        kind===:odd && (n*=z)
        peval(p,x)*d==peval(q,x)*n||error("independent exact backward evaluation disagrees")
    end
    true
end

function count_and_bracket(q,side)
    u=real.(q)
    if side=="even_squared"
        return count_unit(u[1:2:end])
    elseif side=="negative"
        u=u.*[(-1)^(j-1) for j in eachindex(u)]
    end
    count_unit(u)
end

function certify_bound(p,q,side,τ; sine=false)
    s=side=="negative" ? -1 : 1
    pp=real.(p).*[s^(j-1) for j in eachindex(p)]
    qq=real.(q).*[s^(j-1) for j in eachindex(q)]
    # Each half is parameterized by x∈[0,1], with target x or sin(s*x).
    target=sine ? s.*SINE_TAYLOR : QQ[0,1]
    budget=sine ? τ-SINE_TAYLOR_REMAINDER : τ
    err=padd(pp,-pmul(qq,target))
    [(;side,threshold=string(τ),sign,strict_positive(padd(budget.*qq,sign.*err))...) for sign in (-1,1)]
end

is_even_polynomial(p)=all(j->iseven(j-1)||iszero(p[j]),eachindex(p))
function certify_model(rows,axis,uniform,brackets,key,entry,archive_sha)
    case,grid,method=string.(key)
    kind=entry.kind
    kind in (:QTCF,:forced,:native,:even,:odd)||
        error("Unsupported model kind: $kind")
    r=entry.r
    m=exact_model(r,kind)
    phase=case in ("real","odd_sine") ? one(EQ) : case=="imaginary" ? EQ(im) : error("Unrecognized case: $case")
    sine=case=="odd_sine"
    p,q,nonneg=axis_polynomials(m,phase)
    check_model(r,kind,p,q,phase)
    both_real=all(isreal,p)&&all(isreal,q)
    all_even=is_even_polynomial(p)&&is_even_polynomial(q)
    info=(;case,grid,method,kind=string(kind),archive_sha,
        physical_nominal_nodes=hasproperty(entry,:physical_nodes) ? entry.physical_nodes : -1,
        stored_nodes=length(nodes(r)),degree_p=degree(p),degree_q=degree(q),degree_gcd=degree(m.g),
        exact_backward_checks=4,real_after_rotation=both_real,even_after_rotation=all_even,
        q_nonnegative_even_coefficients=nonneg)
    if !both_real
        push!(rows,(;info...,modular_coprime_prime=0,axis_poles=-1,
            axis_status="unresolved_complex_axis_coefficients",uniform_bound=""))
        return
    end
    # Reduce via the determinant product before this separate modular check.
    cp=modular_coprime(p,q)
    polecount=0;counts_done=true
    count_sides=is_even_polynomial(q) ? ("even_squared",) : ("positive","negative")
    for side in count_sides
        a=count_and_bracket(q,side)
        counts_done &= a.rootcount>=0
        n=a.rootcount
        physical=n<0 ? -1 : side=="even_squared" ? 2n-Int(QQ(0) in a.endpoints) : n
        polecount+=physical
        push!(axis,(;case,grid,method,side,exact_root_count=n,physical_axis_poles=physical,
            status=a.status,initial_variations=a.initial_variations,
            visited_subintervals=a.nodes,subdivision_depth=a.maxdepth,
            exact_endpoint_roots=length(a.endpoints)))
        # Preserve the exact isolating dyadic brackets. A variation-one
        # bracket proves one simple root in its open interval.
        for (lo,hi) in vcat(a.brackets,[(x,x) for x in sort(a.endpoints)])
            push!(brackets,(;case,grid,method,side,coordinate=side=="even_squared" ? "x^2" : "|x|",
                left=string(lo),right=string(hi),exact_endpoint=lo==hi))
        end
    end
    if length(count_sides)==2 && iszero(first(q));polecount-=1;end
    status=counts_done ? "exact_bernstein_descartes" : "unresolved_root_count"
    bound=""
    if counts_done && polecount==0
        @assert real(first(q))>0
        bound_sides=all_even && !sine ? ("positive",) : ("positive","negative")
        thresholds=sine ? SINE_CERT_THRESHOLDS : CERT_THRESHOLDS
        for τ in thresholds
            vals=vcat((certify_bound(p,q,side,τ;sine) for side in bound_sides)...)
            append!(uniform,[(;case,grid,method,a...) for a in vals])
            if all(a->a.status=="strictly_positive",vals)
                bound=string(τ)
                break
            end
            any(a->a.status=="unresolved_budget",vals)&&break
        end
    end
    push!(rows,(;info...,modular_coprime_prime=cp,
        axis_poles=counts_done ? polecount : -1,axis_status=status,uniform_bound=bound))
end

function certify_intervals()
    path=joinpath(DEST,"results.jls")
    bytes=read(path)
    archive_sha=bytes2hex(sha256(bytes))
    fits=deserialize(IOBuffer(bytes))
    rows,axis,uniform,brackets,failures=[NamedTuple[] for _ in 1:5]
    out=joinpath(DEST,"interval_certification")
    mkpath(out)
    for ((slug,method),fit) in sort!(collect(fits);by=first)
        slug in ("real","imaginary","odd_sine") || continue
        if fit.best==0
            push!(failures,(;case=slug,method,error="No completed starting pair and no returned fit"))
            println(stderr,last(failures))
            continue
        end
        key=(slug,"adaptive_continuum",method)
        entry=(;r=fit.r,kind=fit.kind,physical_nodes=length(unique(fit.supports[fit.best])))
        try
            certify_model(rows,axis,uniform,brackets,key,entry,archive_sha)
            println(last(rows))
        catch exception
            push!(failures,(;case=slug,method,error=sprint(showerror,exception)))
            println(stderr,last(failures))
        end
        flush(stdout)
    end
    for (name,table) in (("certificates",rows),("axis_detail",axis),
        ("uniform_detail",uniform),("pole_brackets",brackets),("failures",failures))
        save_cert(out,name*".csv",table)
    end
    println("CONTINUUM INTERVAL CERTIFICATION COMPLETE ",length(rows)," models, ",length(failures)," failures")
end

abspath(PROGRAM_FILE)==(@__FILE__) && certify_intervals()
