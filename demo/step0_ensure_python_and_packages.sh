#!/usr/bin/env bash
set -e

##############################
# Required versions
##############################
REQ_PYTHON_MAJOR=3
REQ_PYTHON_MINOR=8

REQ_NUMPY_MM="1.24"   # major.minor
REQ_PANDAS_MM="2.0"
REQ_SKLEARN_MM="1.4"

REQ_NUMPY_RANGE=">=1.24,<1.25"
REQ_PANDAS_RANGE=">=2.0,<2.1"
REQ_SKLEARN_RANGE=">=1.4,<1.5"


##############################################
# Helper: ask user Y/N
##############################################
ask_yn() {
    local prompt="$1"
    while true; do
        read -p "$prompt (Y/N): " yn
        case $yn in
            [Yy]* ) return 0 ;;
            [Nn]* ) return 1 ;;
            * ) echo "Please answer Y or N." ;;
        esac
    done
}


##############################################
# 1. Detect Python 3.8.x
##############################################
detect_python() {
    echo "=== Checking Python ==="

    if command -v python3.8 >/dev/null 2>&1; then
        PYTHON_BIN="python3.8"
    elif command -v python3 >/dev/null 2>&1; then
        PYTHON_BIN="python3"
    elif command -v python >/dev/null 2>&1; then
        PYTHON_BIN="python"
    else
        echo "Python: NOT found."
        echo "  -> Status: does NOT meet requirement (needs Python 3.8.x)."
        echo
        echo "This script cannot install Python itself in a portable way."
        echo "Please install Python 3.8.x using your OS package manager (apt, yum, brew, pyenv, etc.)"
        exit 1
    fi

    PYVER_FULL=$($PYTHON_BIN -c "import sys; print(sys.version.replace('\n',' '))" 2>/dev/null || echo "unknown")
    PYVER_SHORT=$($PYTHON_BIN -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')" 2>/dev/null || echo "unknown")
    PYMAJOR=${PYVER_SHORT%%.*}
    PYMINOR=${PYVER_SHORT#*.}

    echo "Detected Python: $PYTHON_BIN ($PYVER_FULL)"

    if [[ "$PYMAJOR" -eq "$REQ_PYTHON_MAJOR" && "$PYMINOR" -eq "$REQ_PYTHON_MINOR" ]]; then
        echo "  -> Status: meets requirement (Python ${REQ_PYTHON_MAJOR}.${REQ_PYTHON_MINOR}.x)."
    else
        echo "  -> Status: does NOT meet requirement (needs Python ${REQ_PYTHON_MAJOR}.${REQ_PYTHON_MINOR}.x)."
        if ! ask_yn "Do you still want to continue using this Python ($PYVER_SHORT)?"; then
            echo "Exiting by user request."
            exit 1
        fi
    fi

    echo
}

detect_python


##############################################
# 2. Ensure pip
##############################################
echo "=== Ensuring pip is available ==="
$PYTHON_BIN -m ensurepip --upgrade >/dev/null 2>&1 || true
$PYTHON_BIN -m pip install --upgrade pip >/dev/null 2>&1 || true
echo


##############################################
# 3. Package version checker
##############################################
check_pkg_major_minor() {
    # Args: <module_name> <required_major.minor> <label>
    local module="$1"
    local req_mm="$2"
    local label="$3"

    local ver
    ver=$($PYTHON_BIN -c "import $module; print($module.__version__)" 2>/dev/null || echo "MISSING")

    if [[ "$ver" == "MISSING" ]]; then
        echo "$label: NOT installed."
        echo "  -> Status: does NOT meet requirement (${req_mm}.x)."
        return 1
    fi

    local mm
    mm=$(echo "$ver" | cut -d. -f1-2)

    echo "$label: detected version $ver"

    if [[ "$mm" == "$req_mm" ]]; then
        echo "  -> Status: meets requirement (${req_mm}.x)."
        return 0
    else
        echo "  -> Status: does NOT meet requirement (needs ${req_mm}.x)."
        return 1
    fi
}

install_pkg_if_needed() {
    # Args: <pip_name> <module_name> <required_major.minor> <range_spec> <label>
    local pip_name="$1"
    local module="$2"
    local req_mm="$3"
    local range_spec="$4"
    local label="$5"

    echo "=== Checking $label ==="
    if check_pkg_major_minor "$module" "$req_mm" "$label"; then
        echo
        return 0
    fi

    # If we are here, it's missing or wrong version
    if ask_yn "Install/upgrade $pip_name to version range ${range_spec}?"; then
        echo "Installing $pip_name${range_spec} ..."
        $PYTHON_BIN -m pip install --upgrade "${pip_name}${range_spec}"
        echo
    else
        echo "Skipping $label installation by user request."
        echo
    fi
}


##############################################
# 4. Check & possibly install packages
##############################################
install_pkg_if_needed "numpy"        "numpy"   "$REQ_NUMPY_MM"   "$REQ_NUMPY_RANGE"   "NumPy"
install_pkg_if_needed "pandas"       "pandas"  "$REQ_PANDAS_MM"  "$REQ_PANDAS_RANGE"  "pandas"
install_pkg_if_needed "scikit-learn" "sklearn" "$REQ_SKLEARN_MM" "$REQ_SKLEARN_RANGE" "scikit-learn"

##############################################
# 5. Final summary
##############################################
echo "=== Final environment summary ==="
PY_SUMMARY=$($PYTHON_BIN -c "import sys; print(sys.version.replace('\n',' '))" 2>/dev/null || echo "Not available")
echo "Python:       $PY_SUMMARY"
echo "NumPy:        $($PYTHON_BIN -c 'import numpy; print(numpy.__version__)' 2>/dev/null || echo Not installed)"
echo "pandas:       $($PYTHON_BIN -c 'import pandas; print(pandas.__version__)' 2>/dev/null || echo Not installed)"
echo "scikit-learn: $($PYTHON_BIN -c 'import sklearn; print(sklearn.__version__)' 2>/dev/null || echo Not installed)"
echo "================================="

