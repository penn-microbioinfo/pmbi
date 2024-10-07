{
    description = "A Nix-flake-based development environment for python3.12";

    inputs = {
        nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
        nvimenv.url = "path:/home/amsesk/.config/home-manager/FHSEnv";
    };
    outputs = { self, nixpkgs, nvimenv}:
        let
        system = "x86_64-linux";
    pkgs = nixpkgs;
    in 
    {
        devShells."${system}".default = let
            pkgs = import nixpkgs {
                inherit system;
            };
            env = (pkgs.buildFHSUserEnv {
                name = "pmbi-pipzone";
                targetPkgs = pkgs: (with pkgs; [
                        python312
                        python312Packages.pip
                        python312Packages.virtualenv
                        # cudaPackages.cudatoolkit
                ]);
                runScript = "bash";
                }).env; 
                    in pkgs.mkShell {
            inputsFrom = [ 
                env
                nvimenv.devShells."${system}".default
                ];
        };
    };
}
