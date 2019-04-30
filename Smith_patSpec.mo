package Cardiovascular  
  package Model  
    extends Modelica.Icons.ExamplesPackage;

    package Smith2004  
      extends Modelica.Icons.ExamplesPackage;

      package Parts  
        extends Modelica.Icons.UtilitiesPackage;
        type HydraulicLambda = Real(final quantity = "HydraulicLambda", final unit = "1/m3", displayUnit = "1/ml", nominal = 1e-5, min = 0);
      end Parts;
    end Smith2004;
  end Model;

  package Icons  
    model Runnable_Shallow   annotation(experiment(StopTime = 5, __Dymola_NumberOfIntervals = 5000)); end Runnable_Shallow;
  end Icons;
end Cardiovascular;

package Physiolibrary  "Modelica library for Physiology (version 2.3.2-beta)" 
  extends Modelica.Icons.Package;

  package Hydraulic  "Domain with Pressure and Volumetric Flow" 
    extends Modelica.Icons.Package;

    package Components  
      extends Modelica.Icons.Package;

      model Conductor  "Hydraulic resistor, where conductance=1/resistance" 
        extends Hydraulic.Interfaces.OnePort;
        extends Icons.HydraulicResistor;
        parameter Boolean useConductanceInput = false "=true, if external conductance value is used" annotation(Evaluate = true, HideResult = true);
        parameter Types.HydraulicConductance Conductance = 0 "Hydraulic conductance if useConductanceInput=false";
        Types.RealIO.HydraulicConductanceInput cond(start = Conductance) = c if useConductanceInput;
      protected
        Types.HydraulicConductance c;
      equation
        if not useConductanceInput then
          c = Conductance;
        end if;
        q_in.q = c * (q_in.pressure - q_out.pressure);
      end Conductor;

      model Resistor  
        extends Physiolibrary.Hydraulic.Components.Conductor(final Conductance = 1 / Resistance);
        parameter Physiolibrary.Types.HydraulicResistance Resistance "Hydraulic conductance if useConductanceInput=false";
      end Resistor;

      model ElasticVessel  "Elastic container for blood vessels, bladder, lumens" 
        extends Icons.ElasticBalloon;
        extends SteadyStates.Interfaces.SteadyState(state_start = volume_start);
        Interfaces.HydraulicPort_a q_in;
        parameter Types.Volume volume_start = 1e-11 "Volume start value";
        Types.Volume excessVolume "Additional volume, that generate pressure";
        parameter Boolean useV0Input = false "=true, if zero-pressure-volume input is used" annotation(Evaluate = true, HideResult = true);
        parameter Types.Volume ZeroPressureVolume = 1e-11 "Maximal volume, that does not generate pressure if useV0Input=false";
        parameter Types.Volume CollapsingPressureVolume = 1e-12 "Maximal volume, which generate negative collapsing pressure";
        Types.RealIO.VolumeInput zeroPressureVolume(start = ZeroPressureVolume) = zpv if useV0Input;
        parameter Boolean useComplianceInput = false "=true, if compliance input is used" annotation(Evaluate = true, HideResult = true);
        parameter Types.HydraulicCompliance Compliance = 1 "Compliance if useComplianceInput=false";
        Types.RealIO.HydraulicComplianceInput compliance(start = Compliance) = c if useComplianceInput;
        parameter Boolean useExternalPressureInput = false "=true, if external pressure input is used" annotation(Evaluate = true, HideResult = true);
        parameter Types.Pressure ExternalPressure = 0 "External pressure. Set zero if internal pressure is relative to external. Valid only if useExternalPressureInput=false.";
        parameter Types.Pressure MinimalCollapsingPressure = -101325;
        Types.RealIO.PressureInput externalPressure(start = ExternalPressure) = ep if useExternalPressureInput;
        Types.RealIO.VolumeOutput volume;
      protected
        Types.Volume zpv;
        Types.HydraulicCompliance c;
        Types.Pressure ep;
        parameter Types.Pressure a = MinimalCollapsingPressure / log(Modelica.Constants.eps);
      equation
        if not useV0Input then
          zpv = ZeroPressureVolume;
        end if;
        if not useComplianceInput then
          c = Compliance;
        end if;
        if not useExternalPressureInput then
          ep = ExternalPressure;
        end if;
        excessVolume = max(0, volume - zpv);
        q_in.pressure = smooth(0, if noEvent(volume > CollapsingPressureVolume) then excessVolume / c + ep else a * log(max(Modelica.Constants.eps, volume / CollapsingPressureVolume)) + ep);
        state = volume;
        change = q_in.q;
      end ElasticVessel;

      model ElasticVesselElastance  
        extends Physiolibrary.Hydraulic.Components.ElasticVessel(final Compliance = 1 / Elastance);
        parameter Physiolibrary.Types.HydraulicElastance Elastance = 1 "Elastance if useComplianceInput=false";
      end ElasticVesselElastance;

      model Inertia  "Inertia of the volumetric flow" 
        extends SteadyStates.Interfaces.SteadyState(state_start = volumeFlow_start);
        extends Interfaces.OnePort;
        extends Icons.Inertance;
        parameter Types.VolumeFlowRate volumeFlow_start = 0.3 "Volumetric flow start value";
        parameter Types.HydraulicInertance I "Inertance";
      equation
        state = q_in.q;
        change = (q_in.pressure - q_out.pressure) / I;
      end Inertia;

      model IdealValve  
        extends Interfaces.OnePort;
        Boolean open(start = true) "Switching state";
        Real passableVariable(start = 0, final unit = "1") "Auxiliary variable for actual position on the ideal diode characteristic";
        parameter Types.HydraulicConductance _Gon(final min = 0, displayUnit = "l/(mmHg.min)") = 1.2501026264094e-02 "Forward state-on conductance (open valve conductance)";
        parameter Types.HydraulicConductance _Goff(final min = 0, displayUnit = "l/(mmHg.min)") = 1.2501026264094e-12 "Backward state-off conductance (closed valve conductance)";
        parameter Types.Pressure Pknee(final min = 0) = 0 "Forward threshold pressure";
        parameter Boolean useLimitationInputs = false "=true, if Gon and Goff are from inputs" annotation(Evaluate = true, HideResult = true);
        Types.RealIO.HydraulicConductanceInput Gon(start = _Gon) = gon if useLimitationInputs "open valve conductance = infinity for ideal case";
        Types.RealIO.HydraulicConductanceInput Goff(start = _Goff) = goff if useLimitationInputs "closed valve conductance = zero for ideal case";
      protected
        Types.HydraulicConductance gon;
        Types.HydraulicConductance goff;
        constant Types.Pressure unitPressure = 1;
        constant Types.VolumeFlowRate unitFlow = 1;
      equation
        if not useLimitationInputs then
          gon = _Gon;
          goff = _Goff;
        end if;
        open = passableVariable > Modelica.Constants.eps;
        dp = passableVariable * unitFlow * (if open then 1 / gon else 1) + Pknee;
        volumeFlowRate = passableVariable * unitPressure * (if open then 1 else goff) + goff * Pknee;
      end IdealValve;

      model IdealValveResistance  
        extends Physiolibrary.Hydraulic.Components.IdealValve(final _Gon = 1 / _Ron);
        parameter Physiolibrary.Types.HydraulicResistance _Ron = 79.993432449 "forward state resistance";
      end IdealValveResistance;
    end Components;

    package Sensors  
      extends Modelica.Icons.SensorsPackage;

      model PressureMeasure  "Hydraulic pressure at port" 
        extends Icons.PressureMeasure;
        Interfaces.HydraulicPort_a q_in;
        Types.RealIO.PressureOutput pressure "Pressure";
      equation
        pressure = q_in.pressure;
        q_in.q = 0;
      end PressureMeasure;
    end Sensors;

    package Interfaces  
      extends Modelica.Icons.InterfacesPackage;

      connector HydraulicPort  "Hydraulical connector with pressure and volumetric flow" 
        Types.Pressure pressure "Pressure";
        flow Types.VolumeFlowRate q "Volume flow";
      end HydraulicPort;

      connector HydraulicPort_a  "Hydraulical inflow connector" 
        extends HydraulicPort;
      end HydraulicPort_a;

      connector HydraulicPort_b  "Hydraulical outflow connector" 
        extends HydraulicPort;
      end HydraulicPort_b;

      partial model OnePort  "Hydraulical OnePort" 
        HydraulicPort_a q_in "Volume inflow";
        HydraulicPort_b q_out "Volume outflow";
        Types.VolumeFlowRate volumeFlowRate "Volumetric flow";
        Types.Pressure dp "Pressure gradient";
      equation
        q_in.q + q_out.q = 0;
        volumeFlowRate = q_in.q;
        dp = q_in.pressure - q_out.pressure;
      end OnePort;
    end Interfaces;
  end Hydraulic;

  package SteadyStates  "Dynamic Simulation / Steady State" 
    extends Modelica.Icons.Package;

    package Interfaces  
      extends Modelica.Icons.InterfacesPackage;

      partial model SteadyState  "Abstract class for any dynamic state calculation (for any derivation), which is driven by SimulationType option." 
        parameter Types.SimulationType Simulation = Types.SimulationType.NormalInit "Dynamic with Initialization or Steady State" annotation(Evaluate = true, HideResult = true);
        parameter Real state_start "State start or init value" annotation(HideResult = true);
        Real state(start = state_start, stateSelect = StateSelect.prefer) "This state must be connected in inherited class definition" annotation(HideResult = true);
        Real change "Dynamic change of state value per minute" annotation(HideResult = true);
      initial equation
        if Simulation == Types.SimulationType.NormalInit then
          state = state_start;
        end if;
      equation
        der(state) = change;
      end SteadyState;
    end Interfaces;
  end SteadyStates;

  package Icons  "Icons for physiological models" 
    extends Modelica.Icons.Package;

    class ElasticBalloon  end ElasticBalloon;

    partial class HydraulicResistor  end HydraulicResistor;

    class PressureMeasure  end PressureMeasure;

    class Inertance  end Inertance;
  end Icons;

  package Types  "Physiological units with nominals" 
    extends Modelica.Icons.Package;

    package Constants  
      extends Modelica.Icons.SourcesPackage;

      block FrequencyConst  "Constant signal of type Frequency" 
        parameter Types.Frequency k "Constant Frequency output value";
        RealIO.FrequencyOutput y "Frequency constant";
      equation
        y = k;
      end FrequencyConst;

      block PressureConst  "Constant signal of type Pressure" 
        parameter Types.Pressure k "Constant Pressure output value";
        RealIO.PressureOutput y "Pressure constant";
      equation
        y = k;
      end PressureConst;
    end Constants;

    package RealIO  
      extends Modelica.Icons.Package;
      connector PressureInput = input Pressure "input Pressure as connector";
      connector PressureOutput = output Pressure "output Pressure as connector";
      connector VolumeInput = input Volume "input Volume as connector";
      connector VolumeOutput = output Volume "output Volume as connector";
      connector FrequencyInput = input Frequency "input Frequency as connector";
      connector FrequencyOutput = output Frequency "output Frequency as connector";
      connector HydraulicConductanceInput = input HydraulicConductance "input HydraulicConductance as connector";
      connector HydraulicComplianceInput = input HydraulicCompliance "input HydraulicCompliance as connector";
    end RealIO;

    type Time = Modelica.SIunits.Time(displayUnit = "min", nominal = 60);
    type Frequency = Modelica.SIunits.Frequency(displayUnit = "1/min");
    type Mass = Modelica.SIunits.Mass(displayUnit = "g", nominal = 1e-3, min = 0);
    type Pressure = Modelica.SIunits.Pressure(displayUnit = "mmHg", nominal = 133.322387415);
    type Volume = Modelica.SIunits.Volume(displayUnit = "ml", nominal = 1e-6, min = 0);
    type VolumeFlowRate = Modelica.SIunits.VolumeFlowRate(displayUnit = "ml/min", nominal = 1e-6 / 60);
    type HydraulicConductance = Real(final quantity = "HydraulicConductance", final unit = "m3/(Pa.s)", displayUnit = "ml/(mmHg.min)", nominal = 1e-6 / (133.322387415 * 60), min = 0);
    type HydraulicResistance = Real(final quantity = "HydraulicConductance", final unit = "(Pa.s)/m3", displayUnit = "(mmHg.min)/ml", nominal = 1e+6 * 133.322387415 * 60, min = 0);
    type HydraulicCompliance = Real(final quantity = "HydraulicCompliance", final unit = "m3/Pa", displayUnit = "ml/mmHg", nominal = 1e-6 / 133.322387415);
    type HydraulicElastance = Real(final quantity = "HydraulicElastance", final unit = "Pa/m3", displayUnit = "mmHg/ml", nominal = 133.322387415 / 1e-6);
    type HydraulicInertance = Real(final quantity = "HydraulicInertance", final unit = "Pa.s2/m3", displayUnit = "mmHg.min2/ml", nominal = 133.322387415 * 60 ^ 2 / 1e-6);
    type SimulationType = enumeration(NoInit "Use start values only as a guess of state values", NormalInit "Initialization by start values", InitSteadyState "Initialization in Steady State (initial derivations are zeros)", SteadyState "Steady State = Derivations are zeros during simulation") "Initialization or Steady state options (to determine model type before simulating)" annotation(Evaluate = true);
  end Types;
  annotation(version = "2.3.2-beta", versionBuild = 1, versionDate = "2015-09-15", dateModified = "2015-09-15 12:49:00Z"); 
end Physiolibrary;

package FMITest  
  model HemodynamicsSmith_shallow  
    extends Cardiovascular.Icons.Runnable_Shallow;
    .Physiolibrary.Hydraulic.Components.ElasticVesselElastance aorta(ZeroPressureVolume = 0, volume_start = 0.0001241, Elastance = 92165766.41999);
    .Physiolibrary.Hydraulic.Components.ElasticVesselElastance venaCava(ZeroPressureVolume = 0, volume_start = 0.0002952, Elastance(displayUnit = "Pa/m3") = 786602.0857485);
    .Physiolibrary.Hydraulic.Components.IdealValveResistance aorticValve(Pknee = 0, _Ron(displayUnit = "(mmHg.s)/ml") = 2399802.97347);
    .Physiolibrary.Hydraulic.Components.Resistor Rsys(Resistance(displayUnit = "(mmHg.s)/ml") = 145054757.50752);
    .Physiolibrary.Hydraulic.Components.IdealValveResistance tricuspidValve(Pknee = 0, _Ron(displayUnit = "(mmHg.s)/ml") = 3159740.5817355);
    .Physiolibrary.Hydraulic.Components.Inertia Lav(I(displayUnit = "mmHg.s2/ml") = 16250.665802014, volumeFlow_start(displayUnit = "m3/s") = -1.4e-8);
    .Physiolibrary.Hydraulic.Components.Inertia Lpv(I(displayUnit = "mmHg.s2/ml") = 19822.372560862, volumeFlow_start(displayUnit = "m3/s") = -1.9e-9);
    .Physiolibrary.Hydraulic.Components.IdealValveResistance pulmonaryValve(Pknee = 0, _Ron(displayUnit = "(mmHg.s)/ml") = 733273.1307825);
    .Physiolibrary.Hydraulic.Components.ElasticVesselElastance pulmonaryArteries(ZeroPressureVolume = 0, useExternalPressureInput = true, volume_start = 3.904e-05, Elastance(displayUnit = "Pa/m3") = 49195960.956135);
    .Physiolibrary.Hydraulic.Components.Resistor Rpul(Resistance(displayUnit = "(mmHg.s)/ml") = 20691634.526808);
    .Physiolibrary.Hydraulic.Components.ElasticVesselElastance pulmonaryVeins(ZeroPressureVolume = 0, useExternalPressureInput = true, volume_start = 0.0008269, Elastance(displayUnit = "Pa/m3") = 973253.4281295);
    .Physiolibrary.Hydraulic.Components.IdealValveResistance mitralValve(Pknee = 0, _Ron(displayUnit = "(mmHg.s)/ml") = 2106493.721157);
    .Physiolibrary.Hydraulic.Components.Inertia Ltc(I(displayUnit = "mmHg.s2/ml") = 10678.18997523, volumeFlow_start(displayUnit = "m3/s") = 0.0001372);
    .Physiolibrary.Hydraulic.Components.Inertia Lmt(I(displayUnit = "mmHg.s2/ml") = 10261.557514558, volumeFlow_start(displayUnit = "m3/s") = 0.0001141);
    Physiolibrary.Types.Constants.FrequencyConst HR(k = 1.3333333333333);
    Physiolibrary.Types.Constants.PressureConst IntraThoracicPressure(k = -533.28954966);
    VentricularInteraction_flat ventricularInteraction_flat(lambdas(displayUnit = "1/m3") = 435000, lambdarv(displayUnit = "1/m3") = 23000, Essept(displayUnit = "mmHg/ml") = 6499999676.0309, V0peri = 0.0002, Pi0sept = 148.00118226939, Pi0rv = 28.757638965416, Pi0lv = 16.038683206025, Pi0peri = 66.701190423724, Esrv = 77993596.637775, Eslv = 383941811.27772, lambdalv = 33000, lambdaperi = 30000);
  equation
    connect(aorta.q_in, Rsys.q_in);
    connect(Rsys.q_out, venaCava.q_in);
    connect(pulmonaryValve.q_out, pulmonaryArteries.q_in);
    connect(pulmonaryArteries.q_in, Rpul.q_in);
    connect(Rpul.q_out, pulmonaryVeins.q_in);
    connect(aorticValve.q_out, aorta.q_in);
    connect(venaCava.q_in, Ltc.q_in);
    connect(pulmonaryVeins.q_in, Lmt.q_in);
    connect(Lav.q_out, aorticValve.q_in);
    connect(Ltc.q_out, tricuspidValve.q_in);
    connect(Lpv.q_out, pulmonaryValve.q_in);
    connect(mitralValve.q_in, Lmt.q_out);
    connect(tricuspidValve.q_out, ventricularInteraction_flat.rvflow);
    connect(Lpv.q_in, ventricularInteraction_flat.rvflow);
    connect(ventricularInteraction_flat.lvflow, Lav.q_in);
    connect(mitralValve.q_out, Lav.q_in);
    connect(HR.y, ventricularInteraction_flat.HR);
    connect(IntraThoracicPressure.y, ventricularInteraction_flat.Pth);
    connect(pulmonaryArteries.externalPressure, IntraThoracicPressure.y);
    connect(pulmonaryVeins.externalPressure, IntraThoracicPressure.y);
    annotation(experiment(StopTime = 5)); 
  end HemodynamicsSmith_shallow;

  model VentricularInteraction_flat  
    .Physiolibrary.Types.Volume Vsept(start = 0.000002);
    .Physiolibrary.Types.Volume Vperi;
    parameter .Physiolibrary.Types.Volume V0sept = 0.000002;
    parameter .Physiolibrary.Types.Volume V0peri;
    .Physiolibrary.Types.Pressure Psept;
    .Physiolibrary.Types.Pressure Pperi;
    parameter .Physiolibrary.Types.Pressure Pi0sept;
    parameter .Physiolibrary.Types.Pressure Pi0rv;
    parameter .Physiolibrary.Types.Pressure Pi0lv;
    parameter .Physiolibrary.Types.Pressure Pi0peri "peak isovolumic pressure";
    parameter .Physiolibrary.Types.HydraulicElastance Essept;
    parameter .Physiolibrary.Types.HydraulicElastance Esrv;
    parameter .Physiolibrary.Types.HydraulicElastance Eslv "elastance of systole";
    parameter Real A = 1;
    parameter Real B = 80;
    parameter Real CC = 60 / B / 2;
    .Physiolibrary.Types.Time tm;
    discrete .Physiolibrary.Types.Time HP "heart period";
    discrete .Physiolibrary.Types.Time t0 "time of beginning of the cardiac cycle";
    discrete .Physiolibrary.Types.Time ts "duration of systole";
    parameter Cardiovascular.Model.Smith2004.Parts.HydraulicLambda lambdas;
    parameter Cardiovascular.Model.Smith2004.Parts.HydraulicLambda lambdarv;
    parameter Cardiovascular.Model.Smith2004.Parts.HydraulicLambda lambdalv;
    parameter Cardiovascular.Model.Smith2004.Parts.HydraulicLambda lambdaperi;
    Physiolibrary.Hydraulic.Interfaces.HydraulicPort_a rvflow;
    Physiolibrary.Hydraulic.Interfaces.HydraulicPort_a lvflow;
    .Physiolibrary.Types.RealIO.FrequencyInput HR;
    .Physiolibrary.Types.RealIO.PressureInput Pth;
    Physiolibrary.Types.RealIO.VolumeOutput Vlv(start = 0.0001042);
    Physiolibrary.Types.RealIO.VolumeOutput Vrv(start = 0.0001042);
  equation
    tm = time - pre(t0);
    when {initial(), tm >= pre(HP)} then
      HP = 1 / HR;
      t0 = time;
      ts = 0.16 + 0.3 * HP;
    end when;
    Psept = lvflow.pressure - rvflow.pressure;
    Psept = (Vsept - V0sept) * A * exp(-B * (tm - CC) ^ 2) * Essept + (1 - A * exp(-B * (tm - CC) ^ 2)) * Pi0sept * (exp(lambdas * Vsept) - 1);
    rvflow.pressure - Pperi = (Vrv + Vsept) * A * exp(-B * (tm - CC) ^ 2) * Esrv + (1 - A * exp(-B * (tm - CC) ^ 2)) * Pi0rv * (exp(lambdarv * (Vrv + Vsept)) - 1);
    der(Vrv) = rvflow.q;
    lvflow.pressure - Pperi = (Vlv - Vsept) * A * exp(-B * (tm - CC) ^ 2) * Eslv + (1 - A * exp(-B * (tm - CC) ^ 2)) * Pi0lv * (exp(lambdalv * (Vlv - Vsept)) - 1);
    der(Vlv) = lvflow.q;
    Vperi = Vrv + Vlv;
    Pperi = Pth + Pi0peri * (exp(lambdaperi * (Vperi - V0peri)) - 1);
  end VentricularInteraction_flat;

  model Smith_patSpec  
    extends HemodynamicsSmith_shallow(venaCava(volume_start = V_vc0), Ltc(volumeFlow_start = Q_tc0), pulmonaryArteries(volume_start = V_pa0), pulmonaryVeins(volume_start = V_pu0), Lpv(volumeFlow_start(displayUnit = "ml/min") = Q_pv0), Lmt(volumeFlow_start(displayUnit = "ml/min") = Q_mt0), Lav(volumeFlow_start(displayUnit = "ml/min") = Q_av0), aorta(volume_start = V_ao0), ventricularInteraction_flat(Vlv(start = V_lv0), Vrv(start = V_rv0)));
    parameter Physiolibrary.Types.Mass BW = 72.3 "Body weight";
    parameter Modelica.SIunits.Length Hgt = 178 "Body height";
    parameter Boolean isMan = true;
    parameter Physiolibrary.Types.Volume TotBV = if isMan then (0.3669 * (Hgt / 100) ^ 3 + 0.03219 * BW + 0.6041) * 1000 else (0.3561 * (Hgt / 100) ^ 3 + 0.03308 * BW + 0.1833) * 1000 "total blood volume ( Nadler et al. Surgery 51:224,1962)";
    constant Real ml2m3 = 1 / 1000 / 1000;
    constant Real min2s = 1 / 60;
    parameter Physiolibrary.Types.Volume CircBV = 0.30 * TotBV * ml2m3;
    parameter Physiolibrary.Types.Volume V_lv0 = 94.6812 / 1500 * CircBV;
    parameter Physiolibrary.Types.Volume V_rv0 = 90.7302 / 1500 * CircBV;
    parameter Physiolibrary.Types.Volume V_pa0 = 43.0123 / 1500 * CircBV;
    parameter Physiolibrary.Types.Volume V_pu0 = 808.458 / 1500 * CircBV;
    parameter Physiolibrary.Types.Volume V_ao0 = 133.338 / 1500 * CircBV;
    parameter Physiolibrary.Types.Volume V_vc0 = 329.780 / 1500 * CircBV;
    parameter Physiolibrary.Types.VolumeFlowRate Q_mt0 = 245.581 * min2s * ml2m3;
    parameter Physiolibrary.Types.VolumeFlowRate Q_av0 = -1e-3 * min2s * ml2m3;
    parameter Physiolibrary.Types.VolumeFlowRate Q_tc0 = 190.066 * min2s * ml2m3;
    parameter Physiolibrary.Types.VolumeFlowRate Q_pv0 = -1e-3 * min2s * ml2m3;
    Physiolibrary.Hydraulic.Sensors.PressureMeasure pressureMeasure;
    Modelica.Blocks.Interfaces.RealOutput pressure;
  equation
    connect(pressureMeasure.q_in, aorta.q_in);
    connect(pressureMeasure.pressure, pressure);
  end Smith_patSpec;
end FMITest;

package ModelicaServices  "ModelicaServices (OpenModelica implementation) - Models and functions used in the Modelica Standard Library requiring a tool specific implementation" 
  extends Modelica.Icons.Package;

  package Machine  
    extends Modelica.Icons.Package;
    final constant Real eps = 1.e-15 "Biggest number such that 1.0 + eps = 1.0";
    final constant Real small = 1.e-60 "Smallest number such that small and -small are representable on the machine";
    final constant Real inf = 1.e+60 "Biggest Real number such that inf and -inf are representable on the machine";
    final constant Integer Integer_inf = OpenModelica.Internal.Architecture.integerMax() "Biggest Integer number such that Integer_inf and -Integer_inf are representable on the machine";
  end Machine;
  annotation(Protection(access = Access.hide), version = "3.2.2", versionBuild = 0, versionDate = "2016-01-15", dateModified = "2016-01-15 08:44:41Z"); 
end ModelicaServices;

package Modelica  "Modelica Standard Library - Version 3.2.2" 
  extends Modelica.Icons.Package;

  package Blocks  "Library of basic input/output control blocks (continuous, discrete, logical, table blocks)" 
    extends Modelica.Icons.Package;

    package Interfaces  "Library of connectors and partial models for input/output blocks" 
      extends Modelica.Icons.InterfacesPackage;
      connector RealOutput = output Real "'output Real' as connector";
    end Interfaces;
  end Blocks;

  package Math  "Library of mathematical functions (e.g., sin, cos) and of functions operating on vectors and matrices" 
    extends Modelica.Icons.Package;

    package Icons  "Icons for Math" 
      extends Modelica.Icons.IconsPackage;

      partial function AxisCenter  "Basic icon for mathematical function with y-axis in the center" end AxisCenter;
    end Icons;

    function asin  "Inverse sine (-1 <= u <= 1)" 
      extends Modelica.Math.Icons.AxisCenter;
      input Real u;
      output .Modelica.SIunits.Angle y;
      external "builtin" y = asin(u);
    end asin;

    function exp  "Exponential, base e" 
      extends Modelica.Math.Icons.AxisCenter;
      input Real u;
      output Real y;
      external "builtin" y = exp(u);
    end exp;
  end Math;

  package Constants  "Library of mathematical constants and constants of nature (e.g., pi, eps, R, sigma)" 
    extends Modelica.Icons.Package;
    final constant Real pi = 2 * Math.asin(1.0);
    final constant Real eps = ModelicaServices.Machine.eps "Biggest number such that 1.0 + eps = 1.0";
    final constant .Modelica.SIunits.Velocity c = 299792458 "Speed of light in vacuum";
    final constant Real mue_0(final unit = "N/A2") = 4 * pi * 1.e-7 "Magnetic constant";
  end Constants;

  package Icons  "Library of icons" 
    extends Icons.Package;

    partial package ExamplesPackage  "Icon for packages containing runnable examples" 
      extends Modelica.Icons.Package;
    end ExamplesPackage;

    partial package Package  "Icon for standard packages" end Package;

    partial package InterfacesPackage  "Icon for packages containing interfaces" 
      extends Modelica.Icons.Package;
    end InterfacesPackage;

    partial package SourcesPackage  "Icon for packages containing sources" 
      extends Modelica.Icons.Package;
    end SourcesPackage;

    partial package SensorsPackage  "Icon for packages containing sensors" 
      extends Modelica.Icons.Package;
    end SensorsPackage;

    partial package UtilitiesPackage  "Icon for utility packages" 
      extends Modelica.Icons.Package;
    end UtilitiesPackage;

    partial package IconsPackage  "Icon for packages containing icons" 
      extends Modelica.Icons.Package;
    end IconsPackage;
  end Icons;

  package SIunits  "Library of type and unit definitions based on SI units according to ISO 31-1992" 
    extends Modelica.Icons.Package;

    package Conversions  "Conversion functions to/from non SI units and type definitions of non SI units" 
      extends Modelica.Icons.Package;

      package NonSIunits  "Type definitions of non SI units" 
        extends Modelica.Icons.Package;
        type Temperature_degC = Real(final quantity = "ThermodynamicTemperature", final unit = "degC") "Absolute temperature in degree Celsius (for relative temperature use SIunits.TemperatureDifference)" annotation(absoluteValue = true);
      end NonSIunits;
    end Conversions;

    type Angle = Real(final quantity = "Angle", final unit = "rad", displayUnit = "deg");
    type Length = Real(final quantity = "Length", final unit = "m");
    type Volume = Real(final quantity = "Volume", final unit = "m3");
    type Time = Real(final quantity = "Time", final unit = "s");
    type Velocity = Real(final quantity = "Velocity", final unit = "m/s");
    type Acceleration = Real(final quantity = "Acceleration", final unit = "m/s2");
    type Frequency = Real(final quantity = "Frequency", final unit = "Hz");
    type Mass = Real(quantity = "Mass", final unit = "kg", min = 0);
    type Pressure = Real(final quantity = "Pressure", final unit = "Pa", displayUnit = "bar");
    type VolumeFlowRate = Real(final quantity = "VolumeFlowRate", final unit = "m3/s");
    type FaradayConstant = Real(final quantity = "FaradayConstant", final unit = "C/mol");
  end SIunits;
  annotation(version = "3.2.2", versionBuild = 3, versionDate = "2016-04-03", dateModified = "2016-04-03 08:44:41Z"); 
end Modelica;

model Smith_patSpec_total
  extends FMITest.Smith_patSpec;
end Smith_patSpec_total;
