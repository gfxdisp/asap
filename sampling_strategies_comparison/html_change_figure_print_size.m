function html_change_figure_print_size( H, width, height )

    oldUnits = get(H,'Units');
    set( H, 'Units', 'centimeters' );

    figPos = get(H,'Position');
    figPos(3) = width;
    figPos(4) = height;
    set(H,'Position', figPos);
    set(H,'PaperPosition', [0 0 width height] );

    set( H, 'Units', oldUnits );

end