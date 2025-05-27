#!/bin/bash

#!/bin/bash

# --- Compile all Fortran sources ---
echo "Compiling Fortran sources in parent directory..."
for f in ../*.f ../*.f90; do
    [ -e "$f" ] || continue  # skip if no such file
    gfortran -c "$f" -o "${f%.f}.o"
    if [ $? -ne 0 ]; then
        echo "Failed to compile $f"
        exit 1
    fi
done

# --- Helper: Map test file to required C++ sources ---
get_cpp_sources() {
    local test_file="$1"
    case "$test_file" in
        *test_LBTConfig.cpp)
            echo "../LBTConfig.cpp"
            ;;
        *test_LBTEMconservation.cpp)
            echo "../LBTConfig.cpp ../main_functions.cpp ../LBTcl.cpp"
            ;;
        # Add more mappings as needed
        *)
            echo ""
            ;;
    esac
}

# --- Compile, run, and clean each test ---
find . -name 'test_*.cpp' | while read TEST_SRC; do
    EXE=$(basename "${TEST_SRC%.cpp}")
    OBJ="${EXE}.o"

    # Compile test source
    echo "Compiling $TEST_SRC ..."
    g++ -std=c++17 -c "$TEST_SRC" -o "$OBJ"
    if [ $? -ne 0 ]; then
        echo "Failed to compile $TEST_SRC"
        continue
    fi

    # Compile required C++ sources
    CPP_SOURCES=$(get_cpp_sources "$TEST_SRC")
    CPP_OBJS=""
    for src in $CPP_SOURCES; do
        OBJNAME=$(basename "${src%.cpp}.o")
        g++ -std=c++17 -c "$src" -o "$OBJNAME"
        if [ $? -ne 0 ]; then
            echo "Failed to compile $src"
            continue 2
        fi
        CPP_OBJS="$CPP_OBJS $OBJNAME"
    done

    # Gather Fortran objects
    F90_OBJS=$(ls ../*.o 2>/dev/null | grep -E '\.o$' || true)

    # Link everything
    echo "Linking $EXE ..."
    g++ -std=c++17 -o "$EXE" "$OBJ" $CPP_OBJS $F90_OBJS -lgfortran
    if [ $? -ne 0 ]; then
        echo "Failed to link $EXE"
        continue
    fi

    # Run test
    echo "Running $EXE ..."
    ./"$EXE"
    RUN_STATUS=$?

    # Clean up
    rm -f "$EXE" "$OBJ" $CPP_OBJS

    if [ $RUN_STATUS -ne 0 ]; then
        echo "Test $EXE failed (exit code $RUN_STATUS)."
    else
        echo "Test $EXE passed."
    fi
    echo "--------------------------------------"
done
