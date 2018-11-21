within ;
package FMITest
  model HemodynamicsSmith_shallow
    extends Cardiovascular.Icons.Runnable_Shallow;
    import Physiolibrary.Hydraulic.Components.*;
    ElasticVesselElastance aorta(
      ZeroPressureVolume=0,
      volume_start=0.0001241,
      Elastance=92165766.41999) annotation (Placement(transformation(extent=
             {{-130,-30},{-110,-10}})));
    ElasticVesselElastance venaCava(
      ZeroPressureVolume=0,
      volume_start=0.0002952,
      Elastance(displayUnit="Pa/m3") = 786602.0857485)
      annotation (Placement(transformation(extent={{-130,24},{-110,44}})));
    IdealValveResistance aorticValve(Pknee=0, _Ron(displayUnit=
            "(mmHg.s)/ml") = 2399802.97347)
      annotation (Placement(transformation(extent={{-68,-30},{-88,-10}})));
    Resistor Rsys(Resistance(displayUnit="(mmHg.s)/ml") = 145054757.50752)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-120,6})));
    IdealValveResistance tricuspidValve(Pknee=0, _Ron(displayUnit=
            "(mmHg.s)/ml") = 3159740.5817355)
      annotation (Placement(transformation(extent={{-62,24},{-42,44}})));
    Inertia Lav(I(displayUnit="mmHg.s2/ml") = 16250.665802014,
        volumeFlow_start(displayUnit="m3/s") = -1.4e-8) annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-44,-20})));
    Inertia Lpv(I(displayUnit="mmHg.s2/ml") = 19822.372560862,
        volumeFlow_start(displayUnit="m3/s") = -1.9e-9)
      annotation (Placement(transformation(extent={{32,24},{52,44}})));
    IdealValveResistance pulmonaryValve(Pknee=0, _Ron(displayUnit=
            "(mmHg.s)/ml") = 733273.1307825)
      annotation (Placement(transformation(extent={{62,24},{82,44}})));
    ElasticVesselElastance pulmonaryArteries(
      ZeroPressureVolume=0,
      useExternalPressureInput=true,
      volume_start=3.904e-05,
      Elastance(displayUnit="Pa/m3") = 49195960.956135)
      annotation (Placement(transformation(extent={{102,24},{122,44}})));
    Resistor Rpul(Resistance(displayUnit="(mmHg.s)/ml") = 20691634.526808)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={114,4})));
    ElasticVesselElastance pulmonaryVeins(
      ZeroPressureVolume=0,
      useExternalPressureInput=true,
      volume_start=0.0008269,
      Elastance(displayUnit="Pa/m3") = 973253.4281295)
      annotation (Placement(transformation(extent={{104,-30},{124,-10}})));
    IdealValveResistance mitralValve(Pknee=0, _Ron(displayUnit=
            "(mmHg.s)/ml") = 2106493.721157)
      annotation (Placement(transformation(extent={{52,-30},{32,-10}})));
    Inertia Ltc(I(displayUnit="mmHg.s2/ml") = 10678.18997523,
        volumeFlow_start(displayUnit="m3/s") = 0.0001372)
      annotation (Placement(transformation(extent={{-88,24},{-68,44}})));
    Inertia Lmt(I(displayUnit="mmHg.s2/ml") = 10261.557514558,
        volumeFlow_start(displayUnit="m3/s") = 0.0001141) annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={70,-20})));
    Physiolibrary.Types.Constants.FrequencyConst HR(k=1.3333333333333)
      annotation (Placement(transformation(extent={{-44,0},{-28,14}})));
    Physiolibrary.Types.Constants.PressureConst IntraThoracicPressure(k=-533.28954966)
      annotation (Placement(transformation(extent={{38,12},{50,20}})));
    VentricularInteraction_flat ventricularInteraction_flat(
      lambdas(displayUnit="1/m3") = 435000,
      lambdarv(displayUnit="1/m3") = 23000,
      Essept(displayUnit="mmHg/ml") = 6499999676.0309,
      V0peri=0.0002,
      Pi0sept=148.00118226939,
      Pi0rv=28.757638965416,
      Pi0lv=16.038683206025,
      Pi0peri=66.701190423724,
      Esrv=77993596.637775,
      Eslv=383941811.27772,
      lambdalv=33000,
      lambdaperi=30000)
      annotation (Placement(transformation(extent={{-18,-12},{20,28}})));

  equation
    connect(aorta.q_in, Rsys.q_in) annotation (Line(
        points={{-120,-20},{-120,-4}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(Rsys.q_out, venaCava.q_in) annotation (Line(
        points={{-120,16},{-120,34}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(pulmonaryValve.q_out, pulmonaryArteries.q_in) annotation (Line(
        points={{82,34},{112,34}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(pulmonaryArteries.q_in, Rpul.q_in) annotation (Line(
        points={{112,34},{114,34},{114,14}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(Rpul.q_out, pulmonaryVeins.q_in) annotation (Line(
        points={{114,-6},{114,-20}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(aorticValve.q_out, aorta.q_in) annotation (Line(
        points={{-88,-20},{-120,-20}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(venaCava.q_in, Ltc.q_in) annotation (Line(
        points={{-120,34},{-88,34}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(pulmonaryVeins.q_in, Lmt.q_in) annotation (Line(
        points={{114,-20},{80,-20}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(Lav.q_out, aorticValve.q_in) annotation (Line(
        points={{-54,-20},{-68,-20}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(Ltc.q_out, tricuspidValve.q_in) annotation (Line(
        points={{-68,34},{-62,34}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(Lpv.q_out, pulmonaryValve.q_in) annotation (Line(
        points={{52,34},{62,34}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(mitralValve.q_in, Lmt.q_out) annotation (Line(
        points={{52,-20},{60,-20}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(tricuspidValve.q_out, ventricularInteraction_flat.rvflow)
      annotation (Line(
        points={{-42,34},{0.62,34},{0.62,28}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(Lpv.q_in, ventricularInteraction_flat.rvflow) annotation (Line(
        points={{32,34},{0.62,34},{0.62,28}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(ventricularInteraction_flat.lvflow, Lav.q_in) annotation (Line(
        points={{1,-12},{2,-12},{2,-20},{-34,-20}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(mitralValve.q_out, Lav.q_in) annotation (Line(
        points={{32,-20},{-34,-20}},
        color={0,0,0},
        thickness=1,
        smooth=Smooth.None));
    connect(HR.y, ventricularInteraction_flat.HR) annotation (Line(
        points={{-26,7},{-22,7},{-22,8},{-14.2,8}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(IntraThoracicPressure.y, ventricularInteraction_flat.Pth)
      annotation (Line(
        points={{51.5,16},{58,16},{58,8},{16.58,8}},
        color={0,190,190},
        smooth=Smooth.None));
    connect(pulmonaryArteries.externalPressure, IntraThoracicPressure.y)
      annotation (Line(
        points={{120,42},{120,46},{86,46},{86,16},{51.5,16}},
        color={0,190,190},
        smooth=Smooth.None));
    connect(pulmonaryVeins.externalPressure, IntraThoracicPressure.y)
      annotation (Line(
        points={{122,-12},{122,16},{51.5,16}},
        color={0,190,190},
        smooth=Smooth.None));
    annotation (
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-160,-100},
              {160,100}}), graphics={Rectangle(extent={{-158,62},{-82,-34}},
            lineColor={28,108,200}),Rectangle(extent={{-84,54},{92,-46}},
            lineColor={28,108,200})}),
      Icon(coordinateSystem(extent={{-160,-100},{160,100}})),
      Documentation(info="<html>
<p>Cardiovascular model implemented per description of Smith et al. Parametrized to the imlpementation of Brian Carlson, UNI MICH</p>
<p>[12] B. W. Smith, J. G. Chase, R. I. Nokes, G. M. Shaw, G. Wake, Minimal Haemodynamic System Model Including Ventricular Interaction and Valve Dynamics., Medical Engineering &amp; Physics 26 (2) (2004) 131&ndash;139. doi:10.1016/j.medengphy.2003.10.001.</p>
<p>[13] CellML implementation at URL: http://models.cellml.org/exposure/9d046663ba5cac5c8a61ac146183614b/smith_chase_nokes_shaw_wake_2004.cellml/view</p>
</html>"),
      experiment(StopTime=5));
  end HemodynamicsSmith_shallow;

  model VentricularInteraction_flat
    import Physiolibrary.Types.*;
    Volume Vsept(start=0.000002), Vperi;
    parameter Volume V0sept=0.000002, V0peri;
    Pressure Psept, Pperi;
    parameter Pressure Pi0sept, Pi0rv, Pi0lv, Pi0peri
      "peak isovolumic pressure";
    parameter HydraulicElastance Essept, Esrv, Eslv
      "elastance of systole";
    parameter Real A=1, B=80, CC=60/B/2;
    Time tm;
    discrete Time HP "heart period";
    discrete Time t0 "time of beginning of the cardiac cycle";
    discrete Time ts "duration of systole";
    parameter Cardiovascular.Model.Smith2004.Parts.HydraulicLambda lambdas;
    parameter Cardiovascular.Model.Smith2004.Parts.HydraulicLambda lambdarv;
    parameter Cardiovascular.Model.Smith2004.Parts.HydraulicLambda lambdalv;
    parameter Cardiovascular.Model.Smith2004.Parts.HydraulicLambda lambdaperi;
    Physiolibrary.Hydraulic.Interfaces.HydraulicPort_a rvflow annotation (
       Placement(transformation(extent={{-48,20},{-28,40}}),
          iconTransformation(extent={{-12,90},{8,110}})));
    Physiolibrary.Hydraulic.Interfaces.HydraulicPort_a lvflow annotation (
       Placement(transformation(extent={{-46,-22},{-26,-2}}),
          iconTransformation(extent={{-10,-110},{10,-90}})));
    RealIO.FrequencyInput HR annotation (Placement(transformation(extent=
              {{-78,-40},{-38,0}}), iconTransformation(extent={{-100,-20},
              {-60,20}})));
    RealIO.PressureInput Pth annotation (Placement(transformation(extent=
              {{-6,24},{14,44}}), iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=180,
          origin={82,0})));
    Physiolibrary.Types.RealIO.VolumeOutput Vlv(start=
          0.0001042) annotation (Placement(
          transformation(extent={{20,-100},{40,-80}}), iconTransformation(extent={
              {20,-100},{40,-80}})));
    Physiolibrary.Types.RealIO.VolumeOutput Vrv(start=0.0001042) annotation (Placement(
          transformation(extent={{20,-100},{40,-80}}), iconTransformation(extent={
              {20,80},{40,100}})));
  equation
    //timing
    tm = time - pre(t0);
    when {initial(),tm >= pre(HP)} then
      HP = 1/HR;
      t0 = time;
      ts = 0.16 + 0.3*HP;
    end when;
    //  septum
    Psept = lvflow.pressure - rvflow.pressure;
    Psept = (Vsept - V0sept)*A*exp(-B*(tm - CC)^2)*Essept + (1 - A*exp(-B
      *(tm - CC)^2))*Pi0sept*(exp(lambdas*Vsept) - 1);
    // rightventricle
    rvflow.pressure - Pperi = (Vrv + Vsept)*A*exp(-B*(tm - CC)^2)*Esrv +
      (1 - A*exp(-B*(tm - CC)^2))*Pi0rv*(exp(lambdarv*(Vrv + Vsept)) - 1);
    der(Vrv) = rvflow.q;
    //leftventricle
    lvflow.pressure - Pperi = (Vlv - Vsept)*A*exp(-B*(tm - CC)^2)*Eslv +
      (1 - A*exp(-B*(tm - CC)^2))*Pi0lv*(exp(lambdalv*(Vlv - Vsept)) - 1);
    der(Vlv) = lvflow.q;
    //pericardium
    Vperi = Vrv + Vlv;
    Pperi = Pth + Pi0peri*(exp(lambdaperi*(Vperi - V0peri)) - 1);
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={
              {-100,-100},{100,100}}), graphics={Text(
                    extent={{102,-32},{76,-20}},
                    lineColor={0,0,255},
                    fillColor={255,170,170},
                    fillPattern=FillPattern.Forward,
                    textString="Pth"),Text(
                    extent={{-100,-22},{-74,-34}},
                    lineColor={0,0,255},
                    textString="HR"),Rectangle(
                    extent={{-20,80},{20,-60}},
                    lineColor={0,0,255},
                    fillPattern=FillPattern.Solid,
                    fillColor={0,0,255}),Text(
                    extent={{-100,-60},{100,-80}},
                    lineColor={0,0,255},
                    textString="%name")}));
  end VentricularInteraction_flat;

  model Smith_patSpec
    extends HemodynamicsSmith_shallow(
      venaCava(volume_start=V_vc0),
      Ltc(volumeFlow_start=Q_tc0),
      pulmonaryArteries(volume_start=V_pa0),
      pulmonaryVeins(volume_start=V_pu0),
      Lpv(volumeFlow_start(displayUnit="ml/min") = Q_pv0),
      Lmt(volumeFlow_start(displayUnit="ml/min") = Q_mt0),
      Lav(volumeFlow_start(displayUnit="ml/min") = Q_av0),
      aorta(volume_start=V_ao0),
      ventricularInteraction_flat(Vlv(start=V_lv0), Vrv(start=V_rv0)));
    parameter Physiolibrary.Types.Mass BW = 72.3 "Body weight";
    parameter Modelica.SIunits.Length Hgt = 178 "Body height";
    parameter Boolean isMan = true;
    parameter Physiolibrary.Types.Volume TotBV = if isMan then ((0.3669 * (Hgt/100)^3) + (0.03219 * BW) + 0.6041) * 1000 else ((0.3561 * (Hgt/100)^3) + (0.03308 * BW) + 0.1833) * 1000 "total blood volume ( Nadler et al. Surgery 51:224,1962)";
    constant Real ml2m3 = 1/1000/1000;
    constant Real min2s = 1/60;
    parameter Physiolibrary.Types.Volume CircBV = 0.30 * TotBV*ml2m3;
    parameter Physiolibrary.Types.Volume V_lv0 = (94.6812/1500) * CircBV;
    parameter Physiolibrary.Types.Volume V_rv0 = (90.7302/1500) * CircBV;
    parameter Physiolibrary.Types.Volume V_pa0 = (43.0123/1500) * CircBV;
    parameter Physiolibrary.Types.Volume V_pu0 = (808.458/1500) * CircBV;
    parameter Physiolibrary.Types.Volume V_ao0 = (133.338/1500) * CircBV;
    parameter Physiolibrary.Types.Volume V_vc0 = (329.780/1500) * CircBV;
    parameter Physiolibrary.Types.VolumeFlowRate Q_mt0 = 245.581*min2s*ml2m3;
    parameter Physiolibrary.Types.VolumeFlowRate Q_av0 = -1e-3*min2s*ml2m3;
    parameter Physiolibrary.Types.VolumeFlowRate Q_tc0 = 190.066*min2s*ml2m3;
    parameter Physiolibrary.Types.VolumeFlowRate Q_pv0 = -1e-3*min2s*ml2m3;
    Physiolibrary.Hydraulic.Sensors.PressureMeasure pressureMeasure
      annotation (Placement(transformation(extent={{-114,-96},{-134,-76}})));
    Modelica.Blocks.Interfaces.RealOutput pressure
      annotation (Placement(transformation(extent={{-136,-100},{-156,-80}})));
  equation
    connect(pressureMeasure.q_in, aorta.q_in) annotation (Line(
        points={{-120,-92},{-120,-20}},
        color={0,0,0},
        thickness=1));
    connect(pressureMeasure.pressure, pressure)
      annotation (Line(points={{-130,-90},{-146,-90}}, color={0,0,127}));
    annotation (Documentation(info="<html>
<p>All params checked manually in the model</p>
<p>Following params are not inlcuded in the Modelica model:</p>
<p><br>&percnt; Left ventricle free wall parameters</p>
<p>V_d_lvf = 0;                                &percnt; LV ES zero P volume (mL)        </p>
<p>V_0_lvf = 0;                                &percnt; LV ED pressure param (mL)    </p>
<p>&percnt; Right ventricle free wall parameters</p>
<p>V_d_rvf = 0;                                &percnt; RV ES zero P volume (mL)</p>
<p>V_0_rvf = 0;                                &percnt; RV ED pressure param (mL)</p>
</html>"));
  end Smith_patSpec;
  annotation (uses(Physiolibrary(version="2.3.2-beta"), Modelica(version=
            "3.2.2")));
end FMITest;
