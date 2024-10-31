# coding: utf-8

import os
import unittest
from os.path import abspath, dirname, join

import numpy as np
import pandas as pd
import sansmic, sansmic.io, sansmic.model, sansmic.libsansmic
import tempfile

testdir = dirname(abspath(str(__file__)))


class TestStageDefinition(unittest.TestCase):
    def test_heights_and_depths(self):
        stage = sansmic.model.StageDefinition(
            brine_injection_height=30,
            brine_production_height=100,
            brine_interface_height=523,
        )
        self.assertEqual(stage.brine_injection_depth, -30)
        self.assertEqual(stage.brine_production_depth, -100)
        self.assertEqual(stage.brine_interface_depth, -523)
        self.assertEqual(stage.blanket_depth, -523)
        stage.brine_interface_depth = 4000 + stage.blanket_depth
        stage.brine_injection_depth = 4000 + stage.brine_injection_depth
        stage.brine_production_depth = 4000 + stage.brine_production_depth
        self.assertEqual(stage.brine_interface_depth, 3477)
        self.assertEqual(stage.brine_production_depth, 3900)
        self.assertEqual(stage.brine_injection_depth, 3970)
        self.assertEqual(stage.blanket_depth, 3477)
        stage.blanket_depth = 3333
        self.assertEqual(stage.blanket_depth, 3333)
        self.assertEqual(stage.brine_interface_depth, 3333)
        with self.assertRaises(ValueError):
            a = stage.brine_production_height + 2
        with self.assertRaises(ValueError):
            a = stage.brine_injection_height + 2
        with self.assertRaises(ValueError):
            a = stage.brine_interface_height + 2
        stage.brine_injection_height = 30
        stage.brine_production_height = -100
        stage.brine_interface_height = 523
        self.assertEqual(stage.brine_injection_depth, -30)
        self.assertEqual(stage.brine_production_depth, -100)
        self.assertEqual(stage.brine_interface_depth, -523)
        stage.brine_injection_height = None
        stage.brine_production_height = None
        stage.brine_interface_height = None
        self.assertIsNone(stage.brine_injection_height)
        self.assertIsNone(stage.brine_production_depth)
        self.assertIsNone(stage.blanket_depth)
        self.assertIsNone(stage.brine_injection_depth)
        self.assertIsNone(stage.brine_production_height)
        self.assertIsNone(stage.brine_interface_height)

    def test_init_options(self):
        stage = sansmic.model.StageDefinition(blanket_depth=3948)
        self.assertEqual(stage.brine_interface_depth, 3948)
        stage = sansmic.model.StageDefinition(
            brine_injection_height=30,
            brine_production_height=100,
            brine_interface_height=523,
        )
        self.assertEqual(stage.brine_injection_height, 30)
        self.assertEqual(stage.brine_production_height, 100)
        self.assertEqual(stage.brine_interface_height, 523)
        with self.assertRaises(TypeError):
            stage = sansmic.model.StageDefinition(defaults=38)
        stage = sansmic.model.StageDefinition(
            defaults={
                "solver_timestep": 0.1,
                "save-frequency": "daily",
            }
        )
        self.assertEqual(stage.save_frequency, "daily")

    def test_custom_setattr(self):
        stage = sansmic.model.StageDefinition()
        stage.simulation_mode = 1
        self.assertEqual(stage.simulation_mode, sansmic.model.SimulationMode.WITHDRAWAL)
        stage.simulation_mode = "leach fill"
        self.assertEqual(stage.simulation_mode, sansmic.model.SimulationMode.LEACH_FILL)
        stage.simulation_mode = sansmic.model.SimulationMode.ORDINARY
        self.assertEqual(stage.simulation_mode, sansmic.model.SimulationMode.ORDINARY)
        stage.set_cavern_sg = -10
        self.assertIsNone(stage.set_cavern_sg)
        stage.stop_condition = 0
        self.assertEqual(stage.stop_condition, sansmic.model.StopCondition.DURATION)
        stage.stop_condition = "volume"
        self.assertEqual(stage.stop_condition, sansmic.model.StopCondition.VOLUME)
        stage.stop_condition = sansmic.model.StopCondition.DEPTH
        self.assertEqual(stage.stop_condition, sansmic.model.StopCondition.DEPTH)
        with self.assertRaises(TypeError):
            stage.stop_condition = 34.28
        with self.assertRaises(TypeError):
            stage.simulation_mode = 14.28
        with self.assertRaises(KeyError):
            stage.stop_condition = "foo"
        with self.assertRaises(KeyError):
            stage.simulation_mode = "bar"
        stage.stop_condition = "  "
        self.assertIsNone(stage.stop_condition)
        stage.simulation_mode = sansmic.libsansmic.CRunMode.WITHDRAWAL
        self.assertEqual(stage.simulation_mode, sansmic.model.SimulationMode.WITHDRAWAL)


class TestEnumFunctions(unittest.TestCase):
    def test_unit_conversion_values(self):
        self.assertAlmostEqual(sansmic.model.Units.inch(), 0.0254, 12)
        self.assertAlmostEqual(sansmic.model.Units.foot(), 0.3048, 12)
        self.assertAlmostEqual(sansmic.model.Units.cubic_foot(), 0.028316846592, 12)
        self.assertAlmostEqual(sansmic.model.Units.barrel(), 0.158987294928, 12)
        self.assertAlmostEqual(sansmic.model.Units.centimeter(), 0.01, 12)
