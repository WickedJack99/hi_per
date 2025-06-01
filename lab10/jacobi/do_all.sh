#!/bin/bash
git pull
make
python3 creator.py
sleep 0.5
squeue
