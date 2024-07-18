within ;
model rocketNozzle
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor C2(C=2472.5)
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-20,-14})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor Ra(R=0.008822)
    annotation (Placement(transformation(extent={{-48,2},{-28,22}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor Rb(R=0.00084755)
    annotation (Placement(transformation(extent={{-12,2},{8,22}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor C4(C=1070.3)
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={54,-14})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor Rd(R=0.098885)
    annotation (Placement(transformation(extent={{64,2},{84,22}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R5_half(R=0.0084759)
    annotation (Placement(transformation(extent={{104,2},{124,22}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R1_half(R=0.0083135)
    annotation (Placement(transformation(extent={{-84,2},{-64,22}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor T1
    annotation (Placement(transformation(extent={{-56,42},{-36,62}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor Rc(R=0.090748)
    annotation (Placement(transformation(extent={{26,2},{46,22}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor T2
    annotation (Placement(transformation(extent={{-20,42},{0,62}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor T3
    annotation (Placement(transformation(extent={{16,42},{36,62}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor T4
    annotation (Placement(transformation(extent={{54,42},{74,62}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor T5
    annotation (Placement(transformation(extent={{92,42},{112,62}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature Ti(T=1273.15)
    annotation (Placement(transformation(extent={{-126,2},{-106,22}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature To(T=293.15)
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={150,12})));
equation
  connect(Ra.port_b, Rb.port_a)
    annotation (Line(points={{-28,12},{-12,12}}, color={191,0,0}));
  connect(C2.port, Rb.port_a)
    annotation (Line(points={{-20,-4},{-20,12},{-12,12}}, color={191,0,0}));
  connect(Rd.port_b, R5_half.port_a)
    annotation (Line(points={{84,12},{104,12}}, color={191,0,0}));
  connect(R1_half.port_b, Ra.port_a)
    annotation (Line(points={{-64,12},{-48,12}}, color={191,0,0}));
  connect(T1.port, Ra.port_a)
    annotation (Line(points={{-56,52},{-56,12},{-48,12}}, color={191,0,0}));
  connect(Rb.port_b, Rc.port_a)
    annotation (Line(points={{8,12},{26,12}}, color={191,0,0}));
  connect(Rc.port_b, Rd.port_a)
    annotation (Line(points={{46,12},{64,12}}, color={191,0,0}));
  connect(T2.port, Rb.port_a)
    annotation (Line(points={{-20,52},{-20,12},{-12,12}}, color={191,0,0}));
  connect(T3.port, Rc.port_a)
    annotation (Line(points={{16,52},{16,12},{26,12}}, color={191,0,0}));
  connect(T4.port, Rd.port_a)
    annotation (Line(points={{54,52},{54,12},{64,12}}, color={191,0,0}));
  connect(C4.port, Rd.port_a)
    annotation (Line(points={{54,-4},{54,12},{64,12}}, color={191,0,0}));
  connect(T5.port, R5_half.port_a)
    annotation (Line(points={{92,52},{92,12},{104,12}}, color={191,0,0}));
  connect(R5_half.port_b, To.port) annotation (Line(points={{124,12},{132,12},{
          132,12},{140,12}}, color={191,0,0}));
  connect(Ti.port, R1_half.port_a)
    annotation (Line(points={{-106,12},{-84,12}}, color={191,0,0}));
  annotation (uses(Modelica(version="4.0.0")), experiment(
      StartTime=1,
      StopTime=60,
      __Dymola_Algorithm="Dassl"));
end rocketNozzle;
