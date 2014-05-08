#!/usr/bin/env sh

tmux has-session -t cagd-unifi-development
if [ $? != 0 ]
then
    tmux new-session -s cagd-unifi-development -n tex -d
    tmux send-keys -t cagd-unifi-development 'emacs cagd.tex -nw' C-m
    # tmux send-keys -t cagd-unifi-development 'vim' C-m
    tmux split-window -h -t cagd-unifi-development
    # tmux select-layout -t cagd-unifi-development main-horizontal
    # tmux send-keys -t cagd-unifi-development:0.1 'make' # here we do not send C-m so we prepare a make command only
    tmux new-window -n bsplines -t cagd-unifi-development
    tmux send-keys -t cagd-unifi-development:1.0 'cd b-splines; emacs bspline.jl -nw' C-m
    tmux split-window -h -t cagd-unifi-development
    tmux send-keys -t cagd-unifi-development:1.1 'cd b-splines; julia' C-m
    tmux split-window -v -t cagd-unifi-development
    tmux send-keys -t cagd-unifi-development:1.2 'cd b-splines' C-m

    tmux select-window -t cagd-unifi-development:1
fi
tmux attach -t cagd-unifi-development
