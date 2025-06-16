using ClosedWaveguideDispersion
using Documenter, DocumenterCitations
using Literate

DocMeta.setdocmeta!(ClosedWaveguideDispersion, :DocTestSetup, :(using ClosedWaveguideDispersion); recursive=true)

# Add bib file
bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

# Generate examples by Literate.jl
function generate_examples(source_dir, output_dir; exclude=[])
    for file in readdir(source_dir)
        if endswith(file, ".jl") && !(file in exclude)
            input = joinpath(source_dir, file)
            output = output_dir
            Literate.markdown(input, output; documenter=true)
        end
    end
end

source = joinpath(@__DIR__, "src", "literate_examples")
output = joinpath(@__DIR__, "src")
generate_examples(source, output)

makedocs(;
    modules=[ClosedWaveguideDispersion],
    authors="jyzhang <peakfind@126.com> and contributors",
    sitename="ClosedWaveguideDispersion.jl",
    format=Documenter.HTML(;
        canonical="https://peakfind.github.io/ClosedWaveguideDispersion.jl",
        edit_link="main",
        assets=String[],
        mathengine=MathJax3(),
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => ["tutorial.md", "homogeneous.md"],
        "Examples" => ["dirichlet.md", "zhang_case2.md", "pa2.md", "pa3.md"],
        "API" => "api.md",
        "References" => "references.md",
    ],
    plugins = [bib],
)

deploydocs(;
    repo="github.com/peakfind/ClosedWaveguideDispersion.jl",
    devbranch="main",
)
