/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/
#ifndef RINGERSELECTORTOOLS_RINGERSELECTORTOOL_H 
#define RINGERSELECTORTOOLS_RINGERSELECTORTOOL_H 



// interfaces
#include "RingerSelectorTools/tools/procedures/IModel.h"
#include "RingerSelectorTools/tools/procedures/IThresholds.h"
#include "RingerSelectorTools/tools/procedures/INormalization.h"

// file header lib
#include "RingerSelectorTools/tools/common/RingerReader.h"

// models lib
#include "RingerSelectorTools/tools/procedures/MultiLayerPerceptron.h"
#include "RingerSelectorTools/tools/procedures/Norm1.h"

#include "AsgTools/AsgMessaging.h"
#include <memory>
#include <string>


namespace Ringer{

  class RingerSelectorTool: public asg::AsgMessaging
  {
  
    public:
      /// Standard constructor
      RingerSelectorTool();
  
      /// Standard destructor
      ~RingerSelectorTool();
  
      StatusCode initialize();
  
      StatusCode finalize();
  
      bool accept( double discriminant, double et, double eta, double mu) const;
  
      double  calculate( double et, double eta, double mu,
                         double deltaeta1, double deltaPoverP, double deltaPhiReescaled,
                         double d0significance, double d0pvunbiased, double eProbabilityHT) const;
  
      
      double  calculate( std::vector<float>& rings, double et, double eta, double mu) const;
  
      double  calculate( std::vector<float>& rings, double et, double eta, double mu,
                         double deltaeta1, double deltaPoverP, double deltaPhiReescaled,
                         double d0significance, double d0pvunbiased, double eProbabilityHT) const;
  
      double  calculate( std::vector<float>& rings, double et, double eta, double mu,
                         double eratio, double reta, double rphi, double rhad, double weta2,
                         double f1, double f3) const;
  
      double  calculate( std::vector<float>& rings, double et, double eta, double mu,
                         double eratio, double reta, double rphi, double rhad, double weta2,
                         double f1, double f3, double deltaeta1, double deltaPoverP,
                         double deltaPhiReescaled, double d0significance, double d0pvunbiased,
                         double eProbabilityHT) const;
  
  
      std::string getWP( void ) const {return m_wp;};
      bool useCaloRings() const {return m_useCaloRings;};
      bool useTrack() const {return m_useTrack;};
      bool useShowerShape() const {return m_useShowerShape;};
      bool useTileCal() const {return m_useTileCal;};
      void setConstantsCalibPath(std::string s){m_constantsCalibPath=s;};
      void setThresholdsCalibPath(std::string s){m_thresholdsCalibPath=s;};
     
      

  

    private:
  
      bool retrieve(double et, double eta, double mu, std::shared_ptr<Ringer::IModel> &discr, 
                                                      std::shared_ptr<Ringer::INormalization> &preproc ) const;
      
      std::vector<std::shared_ptr<Ringer::INormalization>>  m_preprocs;
      std::vector<std::shared_ptr<Ringer::IModel>>          m_discriminators;
      std::vector<std::shared_ptr<Ringer::IThresholds>>     m_cutDefs;
  
      Ringer::RingerReader m_reader;
      
  
      double      m_etCut;
      double      m_lumiCut;
  
      bool m_useTrack;
      bool m_useCaloRings;
      bool m_useShowerShape;
      bool m_removeOutputTansigTF;
      bool m_doPileupCorrection;
      bool m_useTileCal;

      //Discriminator configuration
      std::string m_constantsCalibPath, m_thresholdsCalibPath;
      std::string m_wp; // working point



  
  };// end class

}


#endif
