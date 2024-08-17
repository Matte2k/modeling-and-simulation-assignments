within ;

model Baio232805_Assign2_Part1
  Modelica.Blocks.Sources.Step Propeller_reference(height=210, startTime=5)
    "Angular speed reference for the propeller" annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-190,0})));
  Modelica.Mechanics.Rotational.Components.LossyGear Gear_box(
    useSupport=false,
    ratio=2,
    lossTable=[0,0.99,0.99,0,0; 50,0.98,0.98,0.5,0.5; 100,0.97,0.97,1,1; 210,
        0.96,0.96,1.5,1.5],
    useHeatPort=true)
    annotation (Placement(transformation(extent={{80,-10},{100,10}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Gear_wall_capacitance(C=3000)
              annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={82,-78})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor
    Gear_wall_thermal_conduction(G=100) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={80,-40})));
  Modelica.Mechanics.Rotational.Sources.QuadraticSpeedDependentTorque
    Load_torque(tau_nominal=100, w_nominal=210)
    "Defined using a quadratic speed dependent torque" annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={190,0})));
  Modelica.Mechanics.Rotational.Components.Inertia Shaft_inertia(J=(0.8*0.01)*
        2700*(0.8)^2*(1/12)) "Volume*Rho*L^2*1/12"
    annotation (Placement(transformation(extent={{122,-10},{142,10}})));
  Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={0,60})));
  Modelica.Blocks.Continuous.PID PID(
    k=0.085,
    Ti=0.17,
    Td=0.55)                         annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-110,0})));
  Modelica.Electrical.Analog.Sources.SignalVoltage signalVoltage annotation (
      Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=270,
        origin={-70,0})));
  Modelica.Blocks.Math.Feedback feedback annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-150,0})));
  Modelica.Electrical.Analog.Basic.Resistor Coil_resistance(R=0.1)
    annotation (Placement(transformation(extent={{-60,10},{-40,30}})));
  Modelica.Electrical.Analog.Basic.Inductor Inductance(L=0.01)
    annotation (Placement(transformation(extent={{-20,10},{0,30}})));
  Modelica.Electrical.Analog.Basic.RotationalEMF Motor_constant(k=0.3)
    annotation (Placement(transformation(extent={{0,-10},{20,10}})));
  Modelica.Electrical.Analog.Basic.Ground Ground
    annotation (Placement(transformation(extent={{-40,-40},{-20,-20}})));
  Modelica.Mechanics.Rotational.Components.Inertia Motor_inertia(J=0.001)
    annotation (Placement(transformation(extent={{40,-10},{60,10}})));
equation
  connect(signalVoltage.p, Coil_resistance.p)
    annotation (Line(points={{-70,10},{-70,20},{-60,20}}, color={0,0,255}));
  connect(Coil_resistance.n, Inductance.p)
    annotation (Line(points={{-40,20},{-20,20}}, color={0,0,255}));
  connect(Inductance.n, Motor_constant.p)
    annotation (Line(points={{0,20},{10,20},{10,10}}, color={0,0,255}));
  connect(Motor_constant.flange, Motor_inertia.flange_a)
    annotation (Line(points={{20,0},{40,0}}, color={0,0,0}));
  connect(Motor_inertia.flange_b, Gear_box.flange_a)
    annotation (Line(points={{60,0},{80,0}}, color={0,0,0}));
  connect(Gear_box.heatPort, Gear_wall_thermal_conduction.port_a)
    annotation (Line(points={{80,-10},{80,-30}}, color={191,0,0}));
  connect(Gear_wall_thermal_conduction.port_b, Gear_wall_capacitance.port)
    annotation (Line(points={{80,-50},{80,-68},{82,-68}},
                                                 color={191,0,0}));
  connect(Gear_box.flange_b, Shaft_inertia.flange_a)
    annotation (Line(points={{100,0},{122,0}}, color={0,0,0}));
  connect(Shaft_inertia.flange_b, Load_torque.flange) annotation (Line(points={
          {142,0},{161,0},{161,7.21645e-16},{180,7.21645e-16}}, color={0,0,0}));
  connect(Propeller_reference.y, feedback.u1) annotation (Line(points={{-179,0},
          {-168.5,0},{-168.5,9.99201e-16},{-158,9.99201e-16}}, color={0,0,127}));
  connect(speedSensor.w, feedback.u2)
    annotation (Line(points={{-11,60},{-150,60},{-150,8}}, color={0,0,127}));
  connect(speedSensor.flange, Load_torque.flange) annotation (Line(points={{10,
          60},{160,60},{160,6.66134e-16},{180,6.66134e-16}}, color={0,0,0}));
  connect(feedback.y, PID.u) annotation (Line(points={{-141,-1.16573e-15},{-131,
          -1.16573e-15},{-131,0},{-122,0}}, color={0,0,127}));
  connect(PID.y, signalVoltage.v) annotation (Line(points={{-99,0},{-90,0},{-90,
          2.22045e-15},{-82,2.22045e-15}}, color={0,0,127}));
  connect(Motor_constant.n, signalVoltage.n) annotation (Line(points={{10,-10},
          {10,-20},{-70,-20},{-70,-10}}, color={0,0,255}));
  connect(Ground.p, signalVoltage.n)
    annotation (Line(points={{-30,-20},{-70,-20},{-70,-10}}, color={0,0,255}));
  connect(Motor_constant.n, Ground.p)
    annotation (Line(points={{10,-10},{10,-20},{-30,-20}}, color={0,0,255}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-100},{220,
            100}})),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,-100},{
            220,100}})),
    uses(Modelica(version="4.0.0")),
    experiment(StopTime=120, __Dymola_Algorithm="Radau"));
end Baio232805_Assign2_Part1;
