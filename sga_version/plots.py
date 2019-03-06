
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt

def mix_color(c1, c2, c3=None):
    if c3 is None:
        return (c1[0] + c2[0])/2, (c1[1] + c2[1])/2, (c1[2] + c2[2])/2
    return (c1[0] + c2[0] + c3[0]) / 3, (c1[1] + c2[1] + c3[1]) / 3, (c1[2] + c2[2] + c3[2]) / 3


def plot_venn(sets, names, ax, colors, title=""):
    v = venn3(sets, set_labels=names, ax=ax)
    #c = venn3_circles(sets, ax=ax)
    ax.set_title(title)
    v.get_patch_by_id('100').set_color(colors[names[0]])
    v.get_patch_by_id('010').set_color(colors[names[1]])
    v.get_patch_by_id('001').set_color(colors[names[2]])
    v.get_patch_by_id('100').set_alpha(0.9)
    v.get_patch_by_id('010').set_alpha(0.9)
    v.get_patch_by_id('001').set_alpha(0.9)
    try:
        v.get_patch_by_id('110').set_color(mix_color(colors[names[0]], colors[names[1]]))
        v.get_patch_by_id('011').set_color(mix_color(colors[names[1]], colors[names[2]]))
        v.get_patch_by_id('101').set_color(mix_color(colors[names[0]], colors[names[2]]))
        v.get_patch_by_id('110').set_alpha(1)
        v.get_patch_by_id('011').set_alpha(1)
        v.get_patch_by_id('101').set_alpha(1)
        v.get_patch_by_id('111').set_color(mix_color(colors[names[0]], colors[names[1]], colors[names[2]]))
        v.get_patch_by_id('111').set_alpha(0.9)
    except Exception:
        pass


def plot_histogram(value_lists, labels, colors, minval=0, maxval=None, step=None):
    if not maxval:
        maxval=max(map(max, value_lists))
    plt.rcParams['figure.figsize'] = (15, 6)
    plt.hist([map(lambda v: min(maxval, v), values) for values in value_lists], bins=range(minval, maxval + 1, step),
             density=True, label=labels, color=[colors[x] for x in labels])
    plt.xticks(range(minval, maxval, step))
    plt.legend()




