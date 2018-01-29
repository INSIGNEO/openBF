#!/bin/bash
wget https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.0-linux-x86_64.tar.gz
tar -xzvf julia-0.6.0-linux-x86_64.tar.gz
mv julia-903644385b ~/julia0.6
rm julia-0.6.0-linux-x86_64.tar.gz
echo "alias julia='~/julia0.6/bin/julia'" >> ~/.bashrc
source ~/.bashrc

git clone https://github.com/INSIGNEO/openBF
cd openBF
mkdir -p ~/.julia/v0.6/openBF
cp -r src ~/.julia/v0.6/openBF/
cp main.jl ~/.julia/v0.6/openBF/

mkdir -p ~/.julia/v0.6/BTypes/src
cp src/BTypes.jl ~/.julia/v0.6/BTypes/src/

echo "alias openBF='cp ~/.julia/v0.6/openBF/main.jl ./main.jl && julia main.jl $1'" >> ~/.bashrc
