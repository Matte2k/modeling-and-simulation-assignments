within ;
model Ex2_Part1_simpler
  Modelica.Blocks.Sources.Step Propeller_reference(height=210, startTime=5)
    "Angular speed reference for the propeller" annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={-154,62})));
  Modelica.Mechanics.Rotational.Components.LossyGear Gear_box(
    ratio=2,
    lossTable=[0,0.99,0,0,0; 50,0.98,0,0.5,0; 100,0.97,0,1,0; 210,0.96,0,1.5,0],

    useHeatPort=false)
    annotation (Placement(transformation(extent={{0,-10},{20,10}})));
  Modelica.Mechanics.Rotational.Sources.QuadraticSpeedDependentTorque
    Load_torque(tau_nominal=100, w_nominal=210)
    "Defined using a quadratic speed dependent torque" annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={90,0})));
  Modelica.Mechanics.Rotational.Components.Inertia Shaft_inertia(J=(0.8*0.01)*
        2700*(0.8)^2*(1/12)) "Volume*Rho*L^2*1/12"
    annotation (Placement(transformation(extent={{40,-10},{60,10}})));
  Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={22,30})));
  Modelica.Electrical.Analog.Sources.SignalVoltage signalVoltage annotation (
      Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=270,
        origin={-146,-10})));
  Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-124,-40})));
  Modelica.Mechanics.Rotational.Components.Inertia inertiaRotor(final J=Jr)
    annotation (Placement(transformation(
        origin={-24,0},
        extent={{10,10},{-10,-10}},
        rotation=180)));
  Modelica.Mechanics.Rotational.Components.Fixed fixed if (not useSupport)
    annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=180,
        origin={-72,-46})));
  Modelica.Electrical.Analog.Basic.Resistor ra(final R=0.1, final useHeatPort=
        false)              annotation (Placement(transformation(extent={{10,-10},
            {-10,10}},
        rotation=180,
        origin={-96,40})));
  Modelica.Electrical.Machines.BasicMachines.Components.InductorDC
                                               la(final L=0.01)
    annotation (Placement(transformation(extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-70,40})));
  Modelica.Electrical.Machines.BasicMachines.Components.AirGapDC
                                             airGapDC
                                   annotation (Placement(transformation(extent={{-10,10},
            {10,-10}},           rotation=180,
        origin={-72,-16})));
  Modelica.Electrical.Analog.Basic.Ground eGround annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-100,-36})));
equation
  connect(Load_torque.flange,Shaft_inertia. flange_b)
    annotation (Line(points={{80,0},{60,0}},   color={0,0,0}));
  connect(Shaft_inertia.flange_a,Gear_box. flange_b)
    annotation (Line(points={{40,0},{20,0}},   color={0,0,0}));
  connect(speedSensor.flange,Shaft_inertia. flange_b)
    annotation (Line(points={{32,30},{66,30},{66,0},{60,0}},   color={0,0,0}));
  connect(Propeller_reference.y, signalVoltage.v) annotation (Line(points={{
          -165,62},{-170,62},{-170,-10},{-158,-10}}, color={0,0,127}));
  connect(signalVoltage.n, ground.p) annotation (Line(points={{-146,-20},{-146,
          -40},{-134,-40}}, color={0,0,255}));
  connect(la.p,ra. n)
    annotation (Line(points={{-80,40},{-86,40}},
                                               color={0,0,255}));
  connect(airGapDC.pin_en,eGround. p) annotation (Line(points={{-82,-26},{-82,
          -58},{-110,-58},{-110,-36}},
                                color={0,0,255}));
  connect(airGapDC.pin_ap,la. n) annotation (Line(
      points={{-62,-6},{-58,-6},{-58,26},{-52,26},{-52,40},{-60,40}},
                                color={0,0,255}));
  connect(airGapDC.flange,inertiaRotor. flange_a) annotation (Line(
      points={{-72,-6},{-72,1.72085e-15},{-34,1.72085e-15}}));
  connect(fixed.flange, airGapDC.support)
    annotation (Line(points={{-72,-46},{-72,-26}}, color={0,0,0}));
  connect(Gear_box.flange_a, inertiaRotor.flange_b) annotation (Line(points={{0,
          0},{-7,0},{-7,-7.21645e-16},{-14,-7.21645e-16}}, color={0,0,0}));
  connect(signalVoltage.p, ra.p) annotation (Line(points={{-146,-1.77636e-15},{
          -146,40},{-106,40}}, color={0,0,255}));
  connect(signalVoltage.n, airGapDC.pin_an) annotation (Line(points={{-146,-20},
          {-146,-40},{-140,-40},{-140,-62},{-58,-62},{-58,-26},{-62,-26}},
        color={0,0,255}));
  connect(eGround.p, airGapDC.pin_ep)
    annotation (Line(points={{-110,-36},{-110,-6},{-82,-6}}, color={0,0,255}));
  annotation (
    uses(Modelica(version="4.0.0")),
    Diagram(coordinateSystem(extent={{-160,-100},{160,100}})),
    Icon(coordinateSystem(extent={{-160,-100},{160,100}})));
end Ex2_Part1_simpler;
