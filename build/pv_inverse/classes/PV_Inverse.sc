
PV_DistrN : PV_ChainUGen
{
	*new { arg buffer;
		^this.multiNew('control', buffer)
	}
}


PV_Inverse_0 : PV_ChainUGen
{
	*new { arg buffer, freqC=440.0, freqT=440.0;
		^this.multiNew('control', buffer, freqC, freqT)
	}
}


PV_Inverse_1 : PV_ChainUGen
{
	*new { arg buffer, freqC=440.0, freqT=440.0;
		^this.multiNew('control', buffer, freqC, freqT)
	}
}

PV_Inverse_2 : PV_ChainUGen
{
	*new { arg buffer, freqC=440.0, freqT=440.0;
		^this.multiNew('control', buffer, freqC, freqT)
	}
}


PV_Inverse_Exp : PV_ChainUGen
{
	*new { arg buffer, freqC=440.0, freqT=440.0;
		^this.multiNew('control', buffer, freqC, freqT)
	}
}



PV_Inv : PV_ChainUGen
{
	*new { arg buffer, freqC=440.0 ;
		^this.multiNew('control', buffer, freqC)
	}
}


PV_Transp : PV_ChainUGen
{
	*new { arg buffer, k=1.0 ;
		^this.multiNew('control', buffer, k)
	}
}


PV_InvTransp : PV_ChainUGen
{
	*new { arg buffer, freqC=440.0, freqT=440.0;
		^this.multiNew('control', buffer, freqC, freqT)
	}
}

