#!/bin/bash

(time dftb+)  >& output.dftb+
(time modes) &> output.modes
