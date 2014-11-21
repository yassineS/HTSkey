from cosmos import Execution, Cosmos, rel, Recipe
import tools

def genomkey(execution):
    recipe = Recipe()

    echo = recipe.add_source([tools.Echo(tags={'word': 'hello'}), tools.Echo(tags={'word': 'world'})])
    cat = recipe.add_stage(tools.Cat, echo, rel.One2many([('n', [1, 2])]))
    wc = recipe.add_stage(tools.WordCount, cat)

    execution.run(recipe)


if __name__ == '__main__':
    cosmos_app = Cosmos('sqlite:///sqlite.db', default_drm='lsf')
    cosmos_app.initdb()

    ex = Execution.start(cosmos_app, 'Simple', 'out/simple', max_attempts=3, restart=True, skip_confirm=True)
    genomkey(ex)
