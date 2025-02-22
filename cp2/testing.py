

import argparse

def main():
    
    
    cases = {
        "absorbing": (0.1, 0.9, 0.9),
        "equilibrium": (0.7, 0.5, 0.5),
        "cyclic": (0.8, 0.1, 0.1),
    }
    
    
    parser = argparse.ArgumentParser(description="SIR Model Simulation with CLI Inputs")
    parser.add_argument("-p1", type=float, help="Probability of S -> I transition")
    parser.add_argument("--p2", type=float, help="Probability of I -> R transition")
    parser.add_argument("--p3", type=float, help="Probability of R -> S transition")
    parser.add_argument("--case", choices=["absorbing", "equilibrium", "cyclic"], 
                        help="Choose a predefined case for (p1, p2, p3)")
    
    parser.add_argument("--N", type=int, help="Size of the grid")
    
    args = parser.parse_args()
    
    parser = argparse.ArgumentParser(description="SIR Model Simulation with Different Cases")
    

    print(args.p1, args.p2)
    
main()