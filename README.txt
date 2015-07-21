Cellular Automata with Walsh transform, by Andrew Barrette

Run CAWalsh01.exe to compute/display cellular automata with default parameters. In terminal use format "CAWalsh01.exe [CA rank] [CA width] [interaction radius] [rule]" for more control of the CA size and rule properties (default values are 1 200 2 64).

DISPLAY:
Display window layout is as follows:
top-left: CA space (horizontal) vs time (vertical)
top-right: Spatial (red) and temporal (green) products over 5 pixels
bottom-left: Spatial (red) and temporal (green) averages over 5 pixels
bottom-right: Spatial Walsh transform of top-left data


CONTROLS:
p	Pause
arrows	Left/Right arrows lower/raise the rule
1	Restart CA with delta-function initialization
2	Restart CA with random initialization