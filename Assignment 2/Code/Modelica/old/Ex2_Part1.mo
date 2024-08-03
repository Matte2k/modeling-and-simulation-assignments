within ;
model Ex2_Part1
  Modelica.Blocks.Sources.Step Propeller_reference(height=210, startTime=5)
    "Angular speed reference for the propeller" annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={50,70})));
  Modelica.Mechanics.Rotational.Components.LossyGear Gear_box(ratio=2,
      useHeatPort=true)
    annotation (Placement(transformation(extent={{0,0},{20,20}})));
  Modelica.Electrical.Machines.BasicMachines.DCMachines.DC_PermanentMagnet dcpm(
    Ra=0.1,
    La=0.01,
    Jr=0.001) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-32,0})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Gear_wall_capacitance(
      C=3000) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={0,-50})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor
    Gear_wall_thermal_conduction(G=100) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={0,-20})));
  Modelica.Mechanics.Rotational.Sources.QuadraticSpeedDependentTorque
    Load_torque(tau_nominal=100, w_nominal=210)
    "Defined using a quadratic speed dependent torque" annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={90,10})));
  Modelica.Mechanics.Rotational.Components.Inertia Shaft_inertia(J=(0.8*0.01)*
        2700*(0.8)^2*(1/12)) "Volume*Rho*L^2*1/12"
    annotation (Placement(transformation(extent={{40,0},{60,20}})));
  Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={50,40})));
  Modelica.Blocks.Continuous.PID PID annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-90,30})));
  Modelica.Electrical.Analog.Sources.SignalVoltage signalVoltage annotation (
      Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=270,
        origin={-62,0})));
equation
  connect(Gear_box.heatPort, Gear_wall_thermal_conduction.port_a) annotation (
      Line(points={{0,0},{0,-5},{1.77636e-15,-5},{1.77636e-15,-10}}, color={191,
          0,0}));
  connect(Gear_wall_thermal_conduction.port_b, Gear_wall_capacitance.port)
    annotation (Line(points={{0,-30},{0,-40}}, color={191,0,0}));
  connect(Load_torque.flange, Shaft_inertia.flange_b)
    annotation (Line(points={{80,10},{60,10}}, color={0,0,0}));
  connect(Shaft_inertia.flange_a, Gear_box.flange_b)
    annotation (Line(points={{40,10},{20,10}}, color={0,0,0}));
  connect(speedSensor.flange, Shaft_inertia.flange_b)
    annotation (Line(points={{60,40},{70,40},{70,10},{60,10}}, color={0,0,0}));
  connect(dcpm.flange, Gear_box.flange_a)
    annotation (Line(points={{-32,10},{0,10}}, color={0,0,0}));
  connect(signalVoltage.p, dcpm.pin_ap) annotation (Line(points={{-62,10},{-62,
          12},{-48,12},{-48,6},{-42,6}}, color={0,0,255}));
  connect(signalVoltage.n, dcpm.pin_an) annotation (Line(points={{-62,-10},{-62,
          -12},{-48,-12},{-48,-6},{-42,-6}}, color={0,0,255}));
  connect(PID.y, signalVoltage.v) annotation (Line(points={{-90,19},{-90,
          2.10942e-15},{-74,2.10942e-15}}, color={0,0,127}));
  connect(speedSensor.w, PID.u) annotation (Line(points={{39,40},{-20,40},{-20,
          60},{-90,60},{-90,42}}, color={0,0,127}));
  connect(Propeller_reference.y, PID.u)
    annotation (Line(points={{39,70},{-90,70},{-90,42}}, color={0,0,127}));
  annotation (uses(Modelica(version="4.0.0")));
end Ex2_Part1;
