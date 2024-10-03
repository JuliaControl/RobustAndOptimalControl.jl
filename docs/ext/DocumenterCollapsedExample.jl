"""
Documenter extension module, providing a collapsed example block with code evalulation.
"""
module DocumenterCollapsedExample

import Documenter
using Documenter: MarkdownAST, DOM
using Documenter.HTMLWriter: DCtx
using Documenter.MarkdownAST: Node

# This digs deep into Documenter internals, by defining a new at-block that gets evaluated
# during the Documenter "expansion" step. The expansion of CollapsedExample re-uses the
# standard runner for @example blocks, but creates a custom MarkdownAST block, which then
# is dispatched on in the HTMLWriter (domify).
abstract type CollapsedExample <: Documenter.Expanders.ExpanderPipeline end
Documenter.Selectors.matcher(::Type{CollapsedExample}, node, page, doc) = Documenter.Expanders.iscode(node, "@collapsed-example")
Documenter.Selectors.order(::Type{CollapsedExample}) = 7.9
function Documenter.Selectors.runner(::Type{CollapsedExample}, node, page, doc)
    # The ExampleBlocks runner will fail, unless the code block language is `@example *`,
    # so we override it here.
    node.element.info = "@example"
    Documenter.Selectors.runner(Documenter.Expanders.ExampleBlocks, node, page, doc)
    # Runner will set node.elements to be Documenter.MultiOutput, which we will replace
    # with CollapsedOutput.
    node.element = CollapsedOutput(node.element.codeblock)
    return
end
# This is the MarkdownAST element that replaces Documenter.MultiOutput so that we could
# dispatch on it in the writer.
struct CollapsedOutput <: Documenter.AbstractDocumenterBlock
    codeblock :: MarkdownAST.CodeBlock
end
function Documenter.HTMLWriter.domify(dctx::DCtx, node::Node, ::CollapsedOutput)
    DOM.@tags details summary
    # Documenter.MultiOutput has two types of children nodes: MarkdownAST.CodeBlock
    # is assumed to be the input code block (which gets surrounded by `<details>`),
    # and others should be Documenter.MultiOutputElement, which are the ouputs.
    # We'll use the standard domify() methods for the latter.
    map(node.children) do node
        if node.element isa MarkdownAST.CodeBlock
            details[:style=>"padding: 0rem; border: 0px solid lightgray; color: gray"](
                summary("See code"),
                Documenter.HTMLWriter.domify(dctx, node)
            )
        else
            Documenter.HTMLWriter.domify(dctx, node)
        end
    end
end

end # module
