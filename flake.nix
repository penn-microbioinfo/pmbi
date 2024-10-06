{
  description = "A Nix-flake-based development environment for pmbi";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
  };

  outputs = { self , nixpkgs ,... }: let
    system = "x86_64-linux";
  in {
    devShells."${system}".default = let
      pkgs = import nixpkgs {
        inherit system;
      };
    in pkgs.mkShell {
      packages = with pkgs; [
        neovim
        sqlite
        python3
        nodejs_22
        poetry
      ];

      shellHook = ''
        echo "python `${pkgs.python3}/bin/python --version`"
      '';
    };
  };
}
