using ClosedWaveguideDispersion
using Documenter

DocMeta.setdocmeta!(ClosedWaveguideDispersion, :DocTestSetup, :(using ClosedWaveguideDispersion); recursive=true)

makedocs(;
    modules=[ClosedWaveguideDispersion],
    authors="jyzhang <peakfind@126.com> and contributors",
    sitename="ClosedWaveguideDispersion.jl",
    format=Documenter.HTML(;
        canonical="https://Jiayi Zhang.github.io/ClosedWaveguideDispersion.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Jiayi Zhang/ClosedWaveguideDispersion.jl",
    devbranch="main",
)
