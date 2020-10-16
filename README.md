
# Deterministic Approximation for Submodular Maximization over a Matroid in Nearly Linear Time

This repository is the official implementation of [Deterministic Approximation for Submodular Maximization over a Matroid in Nearly Linear Time](). 

## Application1: Multi-Product Viral Marketing

* **Compilation**

```setup
cd TwinGreedy_MPVM
```

```setup
cmake .
```

```setup
make
```

* **Run** 


> Run the program on the Flixster data set:

```setup
./TwinGreedy -c config_flixster.txt
```

> Run the program on the Epinions data set:

```setup
./TwinGreedy -c config_epinions.txt
```

## Application2: Social Network Monitoring

* **Compilation**

```setup
cd TwinGreedy_SNM
```

```setup
make
```

* **Run**

> Run the program on the Erdos–Rényi random graph

```setup
./er.sh
```

> Run the program on the Barabasi-Albert random graph

```setup
./ba.sh
```

## Contributing

Redistribution and use in source and binary forms, with or without modifications, are permitted for academic purposes, provided that the proper acknowledgements are done.
