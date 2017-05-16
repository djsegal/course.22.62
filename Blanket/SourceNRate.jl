function SourceNRate(P_W)
  # assume P_W is in Watts / m^2
  W_area = 4*pi^2*500*249.5 / 1e4  #wall area in m^2
  P_W *= W_area  #W
  P_W /= 1.602e-19 # eV/s
  SourceNRate = P_W / 14.1e6  #SourceNeutron / s
  SourceNRate
end
