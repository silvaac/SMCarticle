# ACS, 10/2020
# Average abs error function to be used in objective function
function cfun2(a, A)
    cf = mean(abs.(a-A))
    return cf
end
# Permutation function
function  pfun(R)
    Rs = copy(R)
    ts = sample(1:length(R), 2)
    Rs[reverse(ts)] = R[ts]
    return Rs
end
# objective function to be minimized
function efun6(a,b,c,b2,R,nlag=50,n2lag=25)
    A = autocor(R, 1:nlag)
    B = autocor(abs.(R), 1:n2lag)
    C = crosscor(R, abs.(R), -nlag:nlag)
    B2 = autocor(R.*R, 1:n2lag)
    Z = cfun2(a,A)+cfun2(b,B)+cfun2(c,C)+cfun2(b2,B2)
    return Z
end
# main opt loop: adapted from TISEN
function annelOpt(aR,bR,cR,bR2,cmin,rs0,nlag,n2lag,cgoal,maxA)
  temp = 1e-4
  maxT = 20000
  maxG = 2000
  maxS = 200
  wc   = 0.9
  afac = 0.5
  #maxA = 120e6
  rs   = copy(rs0)
  cost = copy(cmin)
  nTot  = 0
  nGood = 0
  nAll  = 0
  ra    = rs
  tini  = temp
  E     = zeros(convert(Int64,ceil(maxA/1e5)),1).+cmin
  tE    = copy(cmin)
  for tt=1:maxA
    nTot = nTot+1
    nAll = nAll+1
    if (nAll%100000) == 0
      println(nAll/1e6)
      E[convert(Int64,ceil(nAll/1e5))] = tE
    end
    cmax = cost-temp*log(rand())
    trs  = pfun(rs)
    tE   = efun6(aR,bR,cR,bR2,trs,nlag,n2lag)
    if tE<=cmax
      nGood = nGood+1
      rs    = trs
      cost  = tE
    end
    if cost<=cgoal
      ra = rs
      break
    end
  if ((nTot<maxT) && (nGood<maxG))
    if (cost<(cmin*wc))
      cmin = cost
      ra   = rs
    end
    continue
  end
  if ((temp==tini) && (nTot>1.5*nGood))
    tini = 10*temp
    temp = tini
  else
    if (nGood<=maxS)
      afac=sqrt(afac)
      maxT=maxT*sqrt(2)
      temp=tini
    else
      temp=temp*afac
    end
  end
  nTot  = 0
  nGood = 0
end
return ra,E
end
