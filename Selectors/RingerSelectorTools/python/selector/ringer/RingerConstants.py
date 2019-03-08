
__all__ = ['Norm1','Sequential']

class SectionClass:
  
  def __init__(self):
    self.EM = 0;
    self.HAD = 1;

Section = SectionClass();

class NormalisationClass:
  
  def __init__(self):
    self.EVENT = 0;
    self.SECTION = 1;
    self.SET = 2;
    self.SEQUENTIAL = 3;
    self.NONE = 4;
    
Normalisation = NormalisationClass();

class LayerClass:
  
  def __init__(self):
    self.PreSamplerB = 0;
    self.EMB1 = 1;
    self.EMB2 = 2;
    self.EMB3 = 3;
    self.PreSamplerE = 4;
    self.EME1 = 5;
    self.EME2 = 6;
    self.EME3 = 7;
    self.HEC0 = 8;
    self.HEC1 = 9;
    self.HEC2 = 10;
    self.HEC3 = 11;
    self.TileBar0 = 12;
    self.TileBar1 = 13;
    self.TileBar2 = 14;
    self.TileGap0 = 15;
    self.TileGap1 = 16;
    self.TileGap2 = 17;
    self.TileExt0 = 18;
    self.TileExt1 = 19;
    self.TileExt2 = 20;
    self.FCAL0 = 21;
    self.FCAL1 = 22;
    self.FCAL2 = 23;
    self.Unknown = 24;

Layer = LayerClass()


class Norm1(object):
  
  def __init__(self):

    # this part config the TrigRingerPrepoc
    self.NormalisationRings = [Normalisation.EVENT,
                               Normalisation.EVENT,
                               Normalisation.EVENT,
                               Normalisation.EVENT,
                               Normalisation.EVENT,
                               Normalisation.EVENT,
                               Normalisation.EVENT];

    self.SectionRings = [Section.EM,
                         Section.EM,
                         Section.EM,
                         Section.EM,
                         Section.HAD,
                         Section.HAD,
                         Section.HAD];

    # Only change this if you know what you are doing
    # default configuration from RingerFex
    self.NRings = [8, 64, 8, 8, 4, 4, 4]

 
class Sequential(object):
  
  def __init__(self):

    # this part config the TrigRingerPrepoc
    self.NormalisationRings = [Normalisation.SEQUENTIAL,
                               Normalisation.SEQUENTIAL,
                               Normalisation.SEQUENTIAL,
                               Normalisation.SEQUENTIAL,
                               Normalisation.SEQUENTIAL,
                               Normalisation.SEQUENTIAL,
                               Normalisation.SEQUENTIAL];

    self.SectionRings = [Section.EM,
                         Section.EM,
                         Section.EM,
                         Section.EM,
                         Section.HAD,
                         Section.HAD,
                         Section.HAD];

    # Only change this if you know what you are doing
    # default configuration from RingerFex
    self.NRings = [8, 64, 8, 8, 4, 4, 4]


