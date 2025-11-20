#!/usr/bin/env bash

echo "NumPy version:     $(python -c 'import numpy as np; print(np.__version__)' 2>/dev/null || echo "Not installed")"
echo "Pandas version:    $(python -c 'import pandas as pd; print(pd.__version__)' 2>/dev/null || echo "Not installed")"
echo "scikit-learn version: $(python -c 'import sklearn; print(sklearn.__version__)' 2>/dev/null || echo "Not installed")"

