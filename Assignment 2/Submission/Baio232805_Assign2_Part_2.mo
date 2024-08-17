within ;
model Baio232805_Assign2_Part2
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow=
        2210.9) "To be set according to part 1"
    annotation (Placement(transformation(extent={{-140,-50},{-120,-30}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor Gear_unit_T
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-90,42})));
  Modelica.Thermal.HeatTransfer.Components.Convection Water_convection
    annotation (Placement(transformation(extent={{-60,-30},{-40,-50}})));
  Modelica.Blocks.Sources.Constant Convective_conductance(k=300) "Water"
    annotation (Placement(transformation(extent={{10,-10},{-10,10}},
        rotation=0,
        origin={-10,-70})));
  Modelica.Thermal.FluidHeatFlow.Components.Pipe pipe(
    medium=Modelica.Thermal.FluidHeatFlow.Media.Water_10degC(),
    m=0,
    T0=278.15,
    T0fixed=false,
    V_flow(start=0, fixed=false),
    frictionLoss=0,
    useHeatPort=true,
    h_g=0)
    "mass of medium is: rho_water * length * pi * radius^2 (999.7*0.4*3.14*0.02^2)"
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={30,-10})));
  Modelica.Thermal.FluidHeatFlow.Components.OpenTank Source_tank(
    medium=Modelica.Thermal.FluidHeatFlow.Media.Water_10degC(),
    T0=278.15,
    T0fixed=true,
    ATank=0.01,
    hTank=0.8,
    pAmbient(displayUnit="Pa"),
    useHeatPort=false,
    level(start=0.8, fixed=true))
    annotation (Placement(transformation(extent={{-60,0},{-40,20}})));
  Modelica.Thermal.FluidHeatFlow.Sources.VolumeFlow volumeFlow(
    medium=Modelica.Thermal.FluidHeatFlow.Media.Water_10degC(),
    m=(999.7*0.4*3.14*0.02^2),
    T0=278.15,
    T0fixed=false,
    V_flow(start=0),                                         useVolumeFlowInput=
       true)
    "mass of medium is: rho_water * length * pi * radius^2 (999.7*0.4*3.14*0.02^2)"
             annotation (Placement(transformation(
        extent={{-20,-20},{0,0}},
        rotation=0)));
  Modelica.Thermal.FluidHeatFlow.Components.OpenTank Storage_tank(
    medium=Modelica.Thermal.FluidHeatFlow.Media.Water_10degC(),
    T0(displayUnit="degC") = 278.15,
    T0fixed=true,
    ATank=0.01,
    hTank=0.8,
    pAmbient(displayUnit="bar") = 100000,
    useHeatPort=true,
    level(start=0.0001, fixed=true))
                                annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={70,10})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
    annotation (Placement(transformation(extent={{80,40},{60,60}})));
  Modelica.Blocks.Logical.Hysteresis hysteresis(
    uLow=273.15 + 5,
    uHigh=273.15 + 10,
    pre_y_start=true)
    annotation (Placement(transformation(extent={{110,-6},{130,14}})));
  Modelica.Blocks.Math.BooleanToReal booleanToReal(realTrue=-6000, realFalse=0)
    annotation (Placement(transformation(extent={{130,40},{110,60}})));
  Modelica.Blocks.Math.BooleanToReal booleanToReal1(realTrue=0.00006, realFalse=
       0) annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-10,50})));
  Modelica.Blocks.Logical.Hysteresis Gear_unit_T_range(uLow=273.15 + 40, uHigh=
        273.15 + 60,
    pre_y_start=true)
                     "keep between 40 and 60"
    annotation (Placement(transformation(extent={{-60,60},{-40,80}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Gear_wall_capacitance(C=3000, T(start=
          376.445, fixed=true))
              annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-90,-72})));
  parameter Modelica.Thermal.FluidHeatFlow.Media.Medium medium=
      Modelica.Thermal.FluidHeatFlow.Media.Water_10degC() "Medium";
equation
  connect(Gear_unit_T_range.y, booleanToReal1.u)
    annotation (Line(points={{-39,70},{-10,70},{-10,62}},
                                                    color={255,0,255}));
  connect(hysteresis.y, booleanToReal.u)
    annotation (Line(points={{131,4},{140,4},{140,50},{132,50}},
                                                   color={255,0,255}));
  connect(fixedHeatFlow.port, Water_convection.solid)
    annotation (Line(points={{-120,-40},{-60,-40}},
                                               color={191,0,0}));
  connect(Water_convection.fluid, pipe.heatPort) annotation (Line(points={{-40,-40},
          {30,-40},{30,-20}},                color={191,0,0}));
  connect(Gear_wall_capacitance.port, Water_convection.solid)
    annotation (Line(points={{-90,-62},{-90,-40},{-60,-40}},
                                                         color={191,0,0}));
  connect(Gear_unit_T.port, Water_convection.solid) annotation (Line(points={{-90,32},
          {-90,-40},{-60,-40}},                color={191,0,0}));
  connect(Convective_conductance.y, Water_convection.Gc)
    annotation (Line(points={{-21,-70},{-50,-70},{-50,-50}},
                                                          color={0,0,127}));
  connect(volumeFlow.flowPort_b, pipe.flowPort_a)
    annotation (Line(points={{0,-10},{20,-10}},color={255,0,0}));
  connect(Source_tank.flowPort, volumeFlow.flowPort_a)
    annotation (Line(points={{-50,0},{-50,-10},{-20,-10}},
                                               color={255,0,0}));
  connect(pipe.flowPort_b, Storage_tank.flowPort)
    annotation (Line(points={{40,-10},{70,-10},{70,0}},   color={255,0,0}));
  connect(Storage_tank.TTank, hysteresis.u)
    annotation (Line(points={{81,4},{108,4}},              color={0,0,127}));
  connect(booleanToReal.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{109,50},
          {80,50}},                                 color={0,0,127}));
  connect(prescribedHeatFlow.port, Storage_tank.heatPort)
    annotation (Line(points={{60,50},{50,50},{50,0},{60,0}},
                                                           color={191,0,0}));
  connect(Gear_unit_T.T, Gear_unit_T_range.u) annotation (Line(points={{-90,53},
          {-90,70},{-62,70}},            color={0,0,127}));
  connect(booleanToReal1.y, volumeFlow.volumeFlow) annotation (Line(points={{-10,39},
          {-10,0}},                               color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-160,-100},{160,
            100}})),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-160,-100},{
            160,100}})),
    uses(Modelica(version="4.0.0")),
    experiment(StopTime=300, __Dymola_Algorithm="Dassl"));
end Baio232805_Assign2_Part2;
