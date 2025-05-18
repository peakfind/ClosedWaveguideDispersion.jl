using ClosedWaveguideDispersion
using Documenter

DocMeta.setdocmeta!(ClosedWaveguideDispersion, :DocTestSetup, :(using ClosedWaveguideDispersion); recursive=true)

makedocs(;
    modules=[ClosedWaveguideDispersion],
    authors="jyzhang <peakfind@126.com> and contributors",
    sitename="ClosedWaveguideDispersion.jl",
    format=Documenter.HTML(;
        canonical="https://peakfind.github.io/ClosedWaveguideDispersion.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => ["tutorial.md", "homogeneous.md"],
        "API" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/peakfind/ClosedWaveguideDispersion.jl",
    devbranch="main",
)
