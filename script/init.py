#!/usr/bin/env python3

# ==============================================================================
# MIT License
#
# Copyright (c) 2019 Princeton University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ==============================================================================


import argparse
import os
import subprocess
from string import Template


def InitProject(project_name):
    """ Initialize the ILA project with the given name """

    # the list of target files (destination path, template file)
    target_list = [
        (os.path.join('cmake', 'FindZ3.cmake'), 'FindZ3.cmake.in'),
        (os.path.join('app', 'main.cc'), 'main.cc.in'),
        (os.path.join('include', project_name,
                      '{}_top.h'.format(project_name)), 'header.h.in'),
        (os.path.join('src', '{}_top.cc'.format(project_name)), 'source.cc.in'),
        ('CMakeLists.txt', 'CMakeLists.txt.top.in')
    ]

    # check if has been initialized
    for target in target_list:
        if os.path.exists(target[0]):
            print('File {} already exists. Stop init.'.format(target[0]))
            return None

    # template fields mapping
    d = {}
    d['project_name'] = project_name
    d['header_guard'] = '{}_TOP_H__'.format(project_name.upper())
    d['header_file_name'] = '{}_top.h'.format(project_name)
    d['source_file_name'] = '{}_top.cc'.format(project_name)
    d['top_method_name'] = 'Get{}Ila'.format(project_name.capitalize())

    # generate
    for target in target_list:
        template_file = os.path.join('template', target[1])
        with open(template_file, 'r') as fr:
            template_raw = fr.read()
            template = Template(template_raw)

        dir_name = os.path.dirname(target[0])
        if dir_name:
            os.makedirs(dir_name)

        with open(target[0], 'w') as fw:
            ctn = template.safe_substitute(d)
            fw.write(ctn)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Initialize the project')
    parser.add_argument('project_name', type=str, help='project name')
    args = parser.parse_args()

    InitProject(args.project_name)
