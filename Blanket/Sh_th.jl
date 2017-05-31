function Sh_th(MagMax, P_W, T)
  #MagMax is from magnet team and is in neutrons / m^2

  a = 7.51e-08
  b = -0.1539

  P_W = 1.06;

  SourceN = SourceNRate(P_W)

  MagMaxFlux = MagMax / (T* 3.154e7) #max neutrons / m^2 /s
  MaxNperSourceN = MagMaxFlux / SourceN / 1e4 #n/SourceN-cm^2

choices = [log(MaxNperSourceN/a) / b, 0]
Sh_th = maximum(choices)

end
