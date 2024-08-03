within ;
model Baio232805_Assign2_Part2
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow=
        2210.9) "To be set according to part 1"
    annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor Gear_unit_T
    annotation (Placement(transformation(extent={{-20,-42},{0,-22}})));
  Modelica.Thermal.HeatTransfer.Components.Convection Water_convection
    annotation (Placement(transformation(extent={{-34,-10},{-14,10}})));
  Modelica.Blocks.Sources.Constant Convective_conductance(k=300) "Water"
    annotation (Placement(transformation(extent={{-60,40},{-40,60}})));
  Modelica.Thermal.FluidHeatFlow.Components.Pipe pipe(
    medium=Modelica.Thermal.FluidHeatFlow.Media.Water_10degC(),
    m=(999.7*0.4*0.02^2),
    T0=283.15,
    T0fixed=false,
    tapT=1,
    V_flow(start=0.00006),
    frictionLoss=0,
    useHeatPort=true,
    h_g=0)
    "mass of medium is: rho_water * length * pi * radius^2 (999.7*0.4*3.14*0.02^2)"
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={60,0})));
  Modelica.Thermal.FluidHeatFlow.Components.OpenTank Source_tank(
    medium=Modelica.Thermal.FluidHeatFlow.Media.Water_10degC(),
    T0=283.15,
    T0fixed=true,
    ATank=0.01,
    hTank=0.8,
    pAmbient(displayUnit="Pa") = 101325,
    useHeatPort=false,
    level(start=0.8, fixed=true))
    annotation (Placement(transformation(extent={{50,60},{70,80}})));
  Modelica.Thermal.FluidHeatFlow.Sources.VolumeFlow volumeFlow(
    medium=Modelica.Thermal.FluidHeatFlow.Media.Water_10degC(),
    m=0,
    T0=283.15,
    T0fixed=false,
    V_flow(start=0.00006),                                   useVolumeFlowInput=
       true) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={60,38})));
  Modelica.Thermal.FluidHeatFlow.Components.OpenTank Storage_tank(
    medium=Modelica.Thermal.FluidHeatFlow.Media.Water_10degC(),
    T0(displayUnit="degC") = 278.15,
    T0fixed=true,
    ATank=0.01,
    hTank=0.8,
    pAmbient(displayUnit="Pa") = 101325,
    useHeatPort=true,
    level(start=0, fixed=true)) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={90,-40})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
    annotation (Placement(transformation(extent={{160,-20},{140,0}})));
  Modelica.Blocks.Logical.Hysteresis hysteresis(
    uLow=273.15 + 5,
    uHigh=273.15 + 10,
    pre_y_start=false)
    annotation (Placement(transformation(extent={{120,-80},{140,-60}})));
  Modelica.Blocks.Math.BooleanToReal booleanToReal(realTrue=-6000, realFalse=0)
    annotation (Placement(transformation(extent={{160,-80},{180,-60}})));
  Modelica.Blocks.Math.BooleanToReal booleanToReal1(realTrue=0.00006, realFalse=
       0) annotation (Placement(transformation(extent={{112,-132},{132,-112}})));
  Modelica.Blocks.Logical.Hysteresis Gear_unit_T_range(uLow=273.15 + 40, uHigh=
        273.15 + 60,
    pre_y_start=false)
                     "keep between 40 and 60"
    annotation (Placement(transformation(extent={{40,-132},{60,-112}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Gear_wall_capacitance(C=3000,
    T(start=376.478, fixed=true),
    der_T(start=0.736967))
              annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-48,-52})));
  parameter Modelica.Thermal.FluidHeatFlow.Media.Medium medium=
      Modelica.Thermal.FluidHeatFlow.Media.Water_10degC() "Medium";
equation
  connect(Gear_unit_T_range.y, booleanToReal1.u)
    annotation (Line(points={{61,-122},{110,-122}}, color={255,0,255}));
  connect(hysteresis.y, booleanToReal.u)
    annotation (Line(points={{141,-70},{158,-70}}, color={255,0,255}));
  connect(fixedHeatFlow.port, Water_convection.solid)
    annotation (Line(points={{-60,0},{-34,0}}, color={191,0,0}));
  connect(Water_convection.fluid, pipe.heatPort) annotation (Line(points={{-14,0},
          {50,0},{50,1.77636e-15}},          color={191,0,0}));
  connect(Gear_wall_capacitance.port, Water_convection.solid)
    annotation (Line(points={{-48,-42},{-48,0},{-34,0}}, color={191,0,0}));
  connect(Gear_unit_T.port, Water_convection.solid) annotation (Line(points={{
          -20,-32},{-48,-32},{-48,0},{-34,0}}, color={191,0,0}));
  connect(Convective_conductance.y, Water_convection.Gc)
    annotation (Line(points={{-39,50},{-24,50},{-24,10}}, color={0,0,127}));
  connect(volumeFlow.flowPort_b, pipe.flowPort_a)
    annotation (Line(points={{60,28},{60,10}}, color={255,0,0}));
  connect(Source_tank.flowPort, volumeFlow.flowPort_a)
    annotation (Line(points={{60,60},{60,48}}, color={255,0,0}));
  connect(pipe.flowPort_b, Storage_tank.flowPort)
    annotation (Line(points={{60,-10},{60,-40},{80,-40}}, color={255,0,0}));
  connect(Storage_tank.TTank, hysteresis.u)
    annotation (Line(points={{84,-51},{84,-70},{118,-70}}, color={0,0,127}));
  connect(booleanToReal.y, prescribedHeatFlow.Q_flow) annotation (Line(points={
          {181,-70},{192,-70},{192,-10},{160,-10}}, color={0,0,127}));
  connect(prescribedHeatFlow.port, Storage_tank.heatPort)
    annotation (Line(points={{140,-10},{80,-10},{80,-30}}, color={191,0,0}));
  connect(Gear_unit_T.T, Gear_unit_T_range.u) annotation (Line(points={{1,-32},
          {34,-32},{34,-122},{38,-122}}, color={0,0,127}));
  connect(booleanToReal1.y, volumeFlow.volumeFlow) annotation (Line(points={{
          133,-122},{218,-122},{218,38},{70,38}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    uses(Modelica(version="4.0.0")));
end Baio232805_Assign2_Part2;
