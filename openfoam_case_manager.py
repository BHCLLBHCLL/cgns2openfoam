#!/usr/bin/env python3
"""
OpenFOAM Case Manager
=====================
A PyQt5-based GUI tool for managing OpenFOAM simulation cases.

Features:
  • Copy OpenFOAM case settings to a new project (excluding polyMesh)
  • Map boundary/initial conditions from source patches to target mesh patches
  • Edit OpenFOAM dictionary files with syntax highlighting
  • Visualize case directory structure

Requirements:
  pip install PyQt5

Usage:
  python3 openfoam_case_manager.py
"""

import sys
import os
import re
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QSplitter, QTreeWidget, QTreeWidgetItem, QPlainTextEdit,
    QLabel, QPushButton, QLineEdit, QFileDialog, QTabWidget, QGroupBox,
    QScrollArea, QMessageBox, QDialog, QDialogButtonBox,
    QTableWidget, QTableWidgetItem, QHeaderView, QToolBar,
    QStatusBar, QComboBox, QAction, QFrame,     QSizePolicy,
    QAbstractItemView, QMenu, QCheckBox, QTextEdit, QProgressDialog,
    QStyleFactory, QFormLayout,
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QSize, QRect, QPoint
from PyQt5.QtGui import (
    QFont, QColor, QSyntaxHighlighter, QTextCharFormat,
    QPainter, QTextCursor, QFontMetrics, QPalette,
)

# ══════════════════════════════════════════════════════════════
#   Syntax Highlighter
# ══════════════════════════════════════════════════════════════

class OpenFOAMHighlighter(QSyntaxHighlighter):
    """Syntax highlighter for OpenFOAM dictionary files."""

    def __init__(self, document):
        super().__init__(document)
        self._rules: List[Tuple] = []

        def fmt(color, bold=False, italic=False):
            f = QTextCharFormat()
            f.setForeground(QColor(color))
            if bold:
                f.setFontWeight(QFont.Bold)
            if italic:
                f.setFontItalic(True)
            return f

        # OpenFOAM dictionary keywords
        keywords = [
            'FoamFile', 'version', 'format', 'class', 'object', 'location',
            'dimensions', 'internalField', 'boundaryField',
            'type', 'value', 'gradient', 'refValue', 'refGradient',
            'valueFraction', 'snGradScheme', 'uniform', 'nonuniform',
            'startTime', 'endTime', 'deltaT', 'writeInterval', 'application',
            'startFrom', 'stopAt', 'writeFormat', 'writePrecision',
            'writeCompression', 'timeFormat', 'timePrecision',
            'runTimeModifiable', 'ddtSchemes', 'gradSchemes', 'divSchemes',
            'laplacianSchemes', 'interpolationSchemes', 'snGradSchemes',
            'fluxRequired', 'solvers', 'PIMPLE', 'PISO', 'SIMPLE',
            'relaxationFactors', 'nNonOrthogonalCorrectors', 'nCorrectors',
            'nOuterCorrectors', 'nCells', 'nFaces', 'startFace',
        ]
        kw_fmt = fmt('#569cd6', bold=True)
        for kw in keywords:
            self._rules.append((re.compile(r'\b' + kw + r'\b'), kw_fmt))

        # BC / solver type values
        bc_types = [
            'fixedValue', 'zeroGradient', 'noSlip', 'slip', 'inletOutlet',
            'outletInlet', 'pressureInletOutletVelocity', 'totalPressure',
            'fixedFluxPressure', 'freestreamPressure', 'freestream',
            'cyclic', 'cyclicAMI', 'symmetryPlane', 'symmetry', 'empty',
            'wedge', 'wall', 'patch', 'mappedWall', 'processor',
            'turbulentIntensityKineticEnergyInlet',
            'turbulentMixingLengthFrequencyInlet',
            'turbulentMixingLengthDissipationRateInlet',
            'epsilonWallFunction', 'omegaWallFunction', 'kqRWallFunction',
            'nutkWallFunction', 'nutWallFunction', 'nutUSpaldingWallFunction',
            'Gauss', 'linear', 'limitedLinear', 'upwind', 'linearUpwind',
            'vanLeer', 'MUSCL', 'corrected', 'uncorrected', 'limited',
            'backward', 'CrankNicolson', 'Euler', 'steadyState',
            'PBiCGStab', 'PCG', 'smoothSolver', 'GAMG',
            'GaussSeidel', 'DIC', 'DILU', 'diagonal', 'ascii', 'binary',
            'simpleFoam', 'pisoFoam', 'pimpleFoam', 'icoFoam', 'sonicFoam',
            'buoyantSimpleFoam', 'buoyantPimpleFoam', 'rhoPimpleFoam',
            'volScalarField', 'volVectorField', 'surfaceScalarField',
            'polyBoundaryMesh', 'pointField', 'labelList',
        ]
        bc_fmt = fmt('#4ec9b0')
        for bt in bc_types:
            self._rules.append((re.compile(r'\b' + bt + r'\b'), bc_fmt))

        # Numbers (integer, float, scientific notation)
        self._rules.append((
            re.compile(r'\b-?\d+\.?\d*([eE][+-]?\d+)?\b'),
            fmt('#b5cea8'),
        ))

        # Quoted strings
        self._rules.append((re.compile(r'"[^"]*"'), fmt('#ce9178')))

        # Braces / brackets
        self._rules.append((re.compile(r'[(){}\[\]]'), fmt('#ffd700')))

        # Semicolons
        self._rules.append((re.compile(r';'), fmt('#808080')))

        # Block comment and line comment formats
        self._block_fmt = fmt('#6a9955', italic=True)
        self._line_fmt = fmt('#6a9955', italic=True)

        # Dimensional brackets [ ... ]
        self._rules.append((re.compile(r'\[[^\]]*\]'), fmt('#c586c0')))

    def highlightBlock(self, text: str):
        # Apply rule-based highlighting
        for pattern, char_fmt in self._rules:
            for m in pattern.finditer(text):
                self.setFormat(m.start(), m.end() - m.start(), char_fmt)

        # Multi-line block comments /* ... */
        self.setCurrentBlockState(0)
        start = 0
        if self.previousBlockState() != 1:
            start = text.find('/*')
        while start >= 0:
            end = text.find('*/', start)
            if end == -1:
                self.setCurrentBlockState(1)
                length = len(text) - start
            else:
                length = end - start + 2
            self.setFormat(start, length, self._block_fmt)
            start = text.find('/*', start + length)

        # Line comments // (only when not inside block comment)
        if self.currentBlockState() != 1:
            m = re.search(r'//.*', text)
            if m:
                self.setFormat(m.start(), len(m.group()), self._line_fmt)


# ══════════════════════════════════════════════════════════════
#   OpenFOAM File Parser
# ══════════════════════════════════════════════════════════════

class OpenFOAMParser:
    """Utilities for reading and modifying OpenFOAM case files."""

    @staticmethod
    def strip_comments(content: str) -> str:
        """Replace comment characters with spaces, preserving positions."""
        chars = list(content)
        i = 0
        while i < len(content):
            if content[i:i+2] == '/*':
                end = content.find('*/', i + 2)
                end = (end + 2) if end != -1 else len(content)
                for j in range(i, end):
                    if chars[j] != '\n':
                        chars[j] = ' '
                i = end
            elif content[i:i+2] == '//':
                end = content.find('\n', i)
                end = end if end != -1 else len(content)
                for j in range(i, end):
                    chars[j] = ' '
                i = end
            else:
                i += 1
        return ''.join(chars)

    @staticmethod
    def find_block_end(content: str, open_pos: int) -> int:
        """Return position of the closing brace matching `content[open_pos]`."""
        depth = 0
        for i in range(open_pos, len(content)):
            if content[i] == '{':
                depth += 1
            elif content[i] == '}':
                depth -= 1
                if depth == 0:
                    return i
        return -1

    @staticmethod
    def _clean_content(raw: str) -> str:
        """
        Return a version of raw suitable for structural parsing.
        Tries comment-stripping first; falls back to raw if key keywords vanish
        (handles files where the block-comment header is improperly closed).
        """
        sc = OpenFOAMParser.strip_comments(raw)
        # Sanity check: if stripped content lost all alpha, fall back
        if re.search(r'[a-zA-Z]', sc):
            return sc
        return raw

    @staticmethod
    def get_boundary_patches(case_path: str) -> List[Dict[str, str]]:
        """
        Read constant/polyMesh/boundary and return list of patch info dicts.
        Each dict contains 'name' and 'type'.
        """
        boundary_file = Path(case_path) / 'constant' / 'polyMesh' / 'boundary'
        if not boundary_file.exists():
            return []
        try:
            raw = boundary_file.read_text(encoding='utf-8', errors='replace')
            sc = OpenFOAMParser._clean_content(raw)

            # Skip past the FoamFile { ... } header
            start = 0
            ff = re.search(r'\bFoamFile\s*\{', sc)
            if ff:
                end = OpenFOAMParser.find_block_end(sc, ff.end() - 1)
                start = (end + 1) if end != -1 else 0

            tail_sc = sc[start:]
            tail_raw = raw[start:]

            # Find the ( ... ) patch list
            p_open = tail_sc.find('(')
            if p_open == -1:
                return []
            depth = 0
            p_close = -1
            for i in range(p_open, len(tail_sc)):
                if tail_sc[i] == '(':
                    depth += 1
                elif tail_sc[i] == ')':
                    depth -= 1
                    if depth == 0:
                        p_close = i
                        break
            if p_close == -1:
                return []

            list_sc = tail_sc[p_open + 1: p_close]
            list_raw = tail_raw[p_open + 1: p_close]

            patches: List[Dict[str, str]] = []
            pos = 0
            while pos < len(list_sc):
                m = re.search(r'\b(\w[\w./-]*)\s*\{', list_sc[pos:])
                if not m:
                    break
                name = m.group(1)
                brace = pos + m.end() - 1
                end = OpenFOAMParser.find_block_end(list_sc, brace)
                if end == -1:
                    break
                block = list_raw[brace: end + 1]
                tm = re.search(r'\btype\s+(\w+)\s*;', block)
                patches.append({
                    'name': name,
                    'type': tm.group(1) if tm else 'patch',
                })
                pos = end + 1
            return patches
        except Exception as exc:
            print(f'Error reading boundary file: {exc}')
            return []

    @staticmethod
    def _parse_patch_blocks(inner_raw: str, inner_sc: str) -> Dict[str, str]:
        """Parse named { } blocks from a boundaryField inner content."""
        blocks: Dict[str, str] = {}
        pos = 0
        while pos < len(inner_sc):
            m = re.search(r'\b(\w[\w./-]*)\s*\{', inner_sc[pos:])
            if not m:
                break
            name = m.group(1)
            brace = pos + m.end() - 1
            if brace >= len(inner_sc) or inner_sc[brace] != '{':
                pos = brace + 1
                continue
            end = OpenFOAMParser.find_block_end(inner_sc, brace)
            if end == -1:
                break
            blocks[name] = inner_raw[brace: end + 1]
            pos = end + 1
        return blocks

    @staticmethod
    def update_boundary_field(
        file_path: str,
        mapping: Dict[str, Optional[str]],
        target_patches: List[Dict[str, str]],
        default_bc_by_type: Dict[str, str],
    ) -> Tuple[bool, str]:
        """
        Rewrite the boundaryField section of an OpenFOAM field file.

        Args:
            file_path:        Path to the field file in the target 0/ directory.
            mapping:          {target_patch_name: source_patch_name or None}
            target_patches:   Ordered list of target patch dicts ('name', 'type').
            default_bc_by_type: {patch_type: bc_type_string} for unmapped patches.
        """
        try:
            raw = Path(file_path).read_text(encoding='utf-8', errors='replace')
            sc = OpenFOAMParser._clean_content(raw)

            bf = re.search(r'\bboundaryField\s*\{', sc)
            if not bf:
                return False, 'No boundaryField found'

            open_pos = bf.end() - 1
            close_pos = OpenFOAMParser.find_block_end(sc, open_pos)
            if close_pos == -1:
                return False, 'Malformed boundaryField block'

            inner_raw = raw[open_pos + 1: close_pos]
            inner_sc = sc[open_pos + 1: close_pos]
            source_blocks = OpenFOAMParser._parse_patch_blocks(inner_raw, inner_sc)

            parts = ['\n']
            for tp in target_patches:
                tname = tp['name']
                ttype = tp['type']
                sname = mapping.get(tname)

                if sname and sname in source_blocks:
                    # Reuse source block, rename patch
                    block = source_blocks[sname]
                    parts.append(f'    {tname}\n    {block}\n')
                else:
                    # Generate a sensible default BC
                    bc_type = default_bc_by_type.get(ttype, 'zeroGradient')
                    if bc_type == 'fixedValue':
                        parts.append(
                            f'    {tname}\n'
                            f'    {{\n'
                            f'        type            fixedValue;\n'
                            f'        value           $internalField;\n'
                            f'    }}\n'
                        )
                    elif bc_type == 'noSlip':
                        parts.append(
                            f'    {tname}\n'
                            f'    {{\n'
                            f'        type            noSlip;\n'
                            f'    }}\n'
                        )
                    else:
                        parts.append(
                            f'    {tname}\n'
                            f'    {{\n'
                            f'        type            {bc_type};\n'
                            f'    }}\n'
                        )

            new_raw = (
                raw[:open_pos + 1]
                + ''.join(parts)
                + raw[close_pos:]
            )
            Path(file_path).write_text(new_raw, encoding='utf-8')
            return True, 'OK'
        except Exception as exc:
            return False, str(exc)

    @staticmethod
    def is_valid_case(path: str) -> bool:
        p = Path(path)
        return p.is_dir() and (p / 'system').is_dir() and (p / 'constant').is_dir()

    @staticmethod
    def get_bc_files(case_path: str) -> Tuple[List[str], str]:
        """Return (list_of_filenames, zero_dir_path) from 0/ or 0.orig/."""
        for zd in ('0', '0.orig'):
            d = Path(case_path) / zd
            if d.is_dir():
                files = sorted(f.name for f in d.iterdir() if f.is_file())
                return files, str(d)
        return [], ''

    @staticmethod
    def get_case_info(case_path: str) -> Dict[str, str]:
        """Return a dict with basic case information."""
        info = {}
        ctrl = Path(case_path) / 'system' / 'controlDict'
        if ctrl.exists():
            raw = ctrl.read_text(encoding='utf-8', errors='replace')
            sc = OpenFOAMParser.strip_comments(raw)
            for key in ('application', 'startTime', 'endTime', 'deltaT'):
                m = re.search(r'\b' + key + r'\s+([^\s;]+)\s*;', sc)
                if m:
                    info[key] = m.group(1)
        patches = OpenFOAMParser.get_boundary_patches(case_path)
        info['patches'] = ', '.join(p['name'] for p in patches) or '(none)'
        info['patch_count'] = str(len(patches))
        return info


# ══════════════════════════════════════════════════════════════
#   Line-Number Gutter + Code Editor
# ══════════════════════════════════════════════════════════════

class _LineNumberArea(QWidget):
    def __init__(self, editor):
        super().__init__(editor)
        self._editor = editor

    def sizeHint(self) -> QSize:
        return QSize(self._editor._line_number_width(), 0)

    def paintEvent(self, event):
        self._editor._paint_line_numbers(event)


class CodeEditor(QPlainTextEdit):
    """QPlainTextEdit with line numbers and OpenFOAM syntax highlighting."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._gutter = _LineNumberArea(self)
        self._hl = OpenFOAMHighlighter(self.document())

        font = QFont('Consolas', 10)
        font.setStyleHint(QFont.Monospace)
        self.setFont(font)
        self.setLineWrapMode(QPlainTextEdit.NoWrap)
        self.setTabStopDistance(QFontMetrics(font).horizontalAdvance(' ') * 4)

        # Dark background
        palette = self.palette()
        palette.setColor(QPalette.Base, QColor('#1e1e1e'))
        palette.setColor(QPalette.Text, QColor('#d4d4d4'))
        self.setPalette(palette)

        self.blockCountChanged.connect(self._update_gutter_width)
        self.updateRequest.connect(self._update_gutter)
        self.cursorPositionChanged.connect(self._highlight_current_line)
        self._update_gutter_width()
        self._highlight_current_line()

    def _line_number_width(self) -> int:
        digits = max(len(str(self.blockCount())), 3)
        return 8 + self.fontMetrics().horizontalAdvance('9') * digits

    def _update_gutter_width(self):
        self.setViewportMargins(self._line_number_width(), 0, 0, 0)

    def _update_gutter(self, rect, dy):
        if dy:
            self._gutter.scroll(0, dy)
        else:
            self._gutter.update(0, rect.y(), self._gutter.width(), rect.height())
        if rect.contains(self.viewport().rect()):
            self._update_gutter_width()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        cr = self.contentsRect()
        self._gutter.setGeometry(QRect(cr.left(), cr.top(), self._line_number_width(), cr.height()))

    def _paint_line_numbers(self, event):
        painter = QPainter(self._gutter)
        painter.fillRect(event.rect(), QColor('#252526'))
        block = self.firstVisibleBlock()
        num = block.blockNumber()
        top = int(self.blockBoundingGeometry(block).translated(self.contentOffset()).top())
        bottom = top + int(self.blockBoundingRect(block).height())
        while block.isValid() and top <= event.rect().bottom():
            if block.isVisible() and bottom >= event.rect().top():
                painter.setPen(QColor('#858585'))
                painter.drawText(
                    0, top, self._gutter.width() - 4,
                    self.fontMetrics().height(),
                    Qt.AlignRight, str(num + 1),
                )
            block = block.next()
            top = bottom
            bottom = top + int(self.blockBoundingRect(block).height())
            num += 1

    def _highlight_current_line(self):
        selections = []
        if not self.isReadOnly():
            sel = QTextEdit.ExtraSelection()
            sel.format.setBackground(QColor('#2d2d30'))
            sel.format.setProperty(QTextCharFormat.FullWidthSelection, True)
            sel.cursor = self.textCursor()
            sel.cursor.clearSelection()
            selections.append(sel)
        self.setExtraSelections(selections)


# ══════════════════════════════════════════════════════════════
#   File Copy Worker (background thread)
# ══════════════════════════════════════════════════════════════

class CopyWorker(QThread):
    """Copies case files from source to target in a background thread."""

    progress = pyqtSignal(str)
    finished = pyqtSignal(bool, str)

    def __init__(self, source: str, target: str):
        super().__init__()
        self.source = Path(source)
        self.target = Path(target)

    def run(self):
        try:
            self.target.mkdir(parents=True, exist_ok=True)
            dirs_to_copy = ['system', '0', '0.orig']
            for d in dirs_to_copy:
                src = self.source / d
                dst = self.target / d
                if src.exists():
                    self.progress.emit(f'Copying {d}/ ...')
                    if dst.exists():
                        shutil.rmtree(str(dst))
                    shutil.copytree(str(src), str(dst))

            # Copy constant/ but skip polyMesh/
            self.progress.emit('Copying constant/ (excluding polyMesh) ...')
            src_const = self.source / 'constant'
            dst_const = self.target / 'constant'
            dst_const.mkdir(exist_ok=True)
            if src_const.exists():
                for item in src_const.iterdir():
                    if item.name.lower() == 'polymesh':
                        self.progress.emit('  Skipping constant/polyMesh/')
                        continue
                    dst_item = dst_const / item.name
                    if item.is_dir():
                        if dst_item.exists():
                            shutil.rmtree(str(dst_item))
                        shutil.copytree(str(item), str(dst_item))
                    else:
                        shutil.copy2(str(item), str(dst_item))

            self.finished.emit(True, 'Case files copied successfully.')
        except Exception as exc:
            self.finished.emit(False, str(exc))


# ══════════════════════════════════════════════════════════════
#   Boundary Condition Mapping Dialog
# ══════════════════════════════════════════════════════════════

NEW_PATCH_LABEL = '(new / default)'
_DEFAULT_BC_OPTIONS = [
    'zeroGradient', 'fixedValue', 'noSlip', 'slip',
    'empty', 'symmetryPlane', 'cyclic', 'wall',
]
_WALL_DEFAULT_BC = {'wall': 'noSlip'}


class BoundaryMappingDialog(QDialog):
    """
    Dialog to map source boundary patches to target boundary patches and
    choose default BC types for unmapped / new target patches.
    """

    def __init__(
        self,
        source_case: str,
        target_case: str,
        source_patches: List[Dict],
        target_patches: List[Dict],
        bc_files: List[str],
        parent=None,
    ):
        super().__init__(parent)
        self.source_case = source_case
        self.target_case = target_case
        self.source_patches = source_patches
        self.target_patches = target_patches
        self.bc_files = bc_files

        self.setWindowTitle('Boundary Condition Mapping')
        self.setMinimumSize(820, 580)
        self._build_ui()
        self._auto_map()

    # ── UI construction ────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setSpacing(10)

        # Header
        header = QLabel(
            '<b>Map source boundary patches to target mesh patches</b><br>'
            'For each <i>target</i> patch select which <i>source</i> patch to '
            'copy boundary conditions from.<br>'
            'Select <tt>(new / default)</tt> to generate a default BC.'
        )
        header.setWordWrap(True)
        root.addWidget(header)

        # Case paths
        paths_box = QGroupBox('Cases')
        paths_layout = QFormLayout(paths_box)
        paths_layout.addRow('Source:', QLabel(self.source_case))
        paths_layout.addRow('Target:', QLabel(self.target_case))
        root.addWidget(paths_box)

        # Patch mapping table
        map_box = QGroupBox('Patch Mapping')
        map_layout = QVBoxLayout(map_box)
        self._table = QTableWidget()
        self._table.setColumnCount(4)
        self._table.setHorizontalHeaderLabels([
            'Target Patch', 'Mesh Type', 'Copy BC from (Source)', 'Default BC Type',
        ])
        self._table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self._table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self._table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)
        self._table.horizontalHeader().setSectionResizeMode(3, QHeaderView.Stretch)
        self._table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._table.setAlternatingRowColors(True)
        self._populate_table()
        map_layout.addWidget(self._table)
        root.addWidget(map_box)

        # BC files selector
        files_box = QGroupBox('Apply to field files (0/)')
        files_layout = QHBoxLayout(files_box)
        self._file_checks: Dict[str, QCheckBox] = {}
        for fname in self.bc_files:
            cb = QCheckBox(fname)
            cb.setChecked(True)
            self._file_checks[fname] = cb
            files_layout.addWidget(cb)
        files_layout.addStretch()
        sel_all = QPushButton('All')
        sel_all.setFixedWidth(50)
        sel_all.clicked.connect(lambda: self._toggle_all_files(True))
        sel_none = QPushButton('None')
        sel_none.setFixedWidth(55)
        sel_none.clicked.connect(lambda: self._toggle_all_files(False))
        files_layout.addWidget(sel_all)
        files_layout.addWidget(sel_none)
        root.addWidget(files_box)

        # Buttons
        btns = QDialogButtonBox()
        apply_btn = btns.addButton('Apply to Selected Files', QDialogButtonBox.AcceptRole)
        btns.addButton(QDialogButtonBox.Cancel)
        apply_btn.setStyleSheet('background:#0e639c; color:white; padding:4px 12px;')
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        root.addWidget(btns)

    def _populate_table(self):
        src_names = [NEW_PATCH_LABEL] + [p['name'] for p in self.source_patches]
        self._table.setRowCount(len(self.target_patches))
        for row, tp in enumerate(self.target_patches):
            # Col 0: target patch name (read-only)
            item = QTableWidgetItem(tp['name'])
            item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            self._table.setItem(row, 0, item)

            # Col 1: mesh type (read-only)
            type_item = QTableWidgetItem(tp['type'])
            type_item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            type_item.setForeground(QColor('#4ec9b0'))
            self._table.setItem(row, 1, type_item)

            # Col 2: source patch dropdown
            src_combo = QComboBox()
            src_combo.addItems(src_names)
            src_combo.currentIndexChanged.connect(
                lambda idx, r=row: self._on_source_changed(r, idx)
            )
            self._table.setCellWidget(row, 2, src_combo)

            # Col 3: default BC type (only active when source == "(new)")
            default_combo = QComboBox()
            default_combo.addItems(_DEFAULT_BC_OPTIONS)
            # Pre-select sensible default
            preset = _WALL_DEFAULT_BC.get(tp['type'], 'zeroGradient')
            idx = default_combo.findText(preset)
            if idx >= 0:
                default_combo.setCurrentIndex(idx)
            self._table.setCellWidget(row, 3, default_combo)

    def _on_source_changed(self, row: int, idx: int):
        """Enable/disable the default BC combo based on source selection."""
        default_combo = self._table.cellWidget(row, 3)
        if default_combo:
            default_combo.setEnabled(idx == 0)  # idx 0 == NEW_PATCH_LABEL

    def _auto_map(self):
        """Auto-map target patches to source patches by exact or partial name match."""
        src_names = [p['name'] for p in self.source_patches]
        for row in range(self._table.rowCount()):
            tname = self._table.item(row, 0).text()
            src_combo = self._table.cellWidget(row, 2)
            if src_combo is None:
                continue
            # Try exact match first, then case-insensitive
            matched = 0
            for i, sname in enumerate(src_names):
                if sname == tname:
                    matched = i + 1  # +1 because index 0 = NEW_PATCH_LABEL
                    break
            if matched == 0:
                for i, sname in enumerate(src_names):
                    if sname.lower() == tname.lower():
                        matched = i + 1
                        break
            src_combo.setCurrentIndex(matched)
            self._on_source_changed(row, matched)

    def _toggle_all_files(self, state: bool):
        for cb in self._file_checks.values():
            cb.setChecked(state)

    # ── Getters ───────────────────────────────────────────────

    def get_mapping(self) -> Dict[str, Optional[str]]:
        """Return {target_patch: source_patch or None}."""
        src_names = [p['name'] for p in self.source_patches]
        mapping: Dict[str, Optional[str]] = {}
        for row in range(self._table.rowCount()):
            tname = self._table.item(row, 0).text()
            src_combo = self._table.cellWidget(row, 2)
            if src_combo is None:
                mapping[tname] = None
                continue
            idx = src_combo.currentIndex()
            mapping[tname] = src_names[idx - 1] if idx > 0 else None
        return mapping

    def get_default_bc_by_type(self) -> Dict[str, str]:
        """Return {patch_type: bc_type} from the unmapped rows."""
        result: Dict[str, str] = {}
        for row in range(self._table.rowCount()):
            ptype = self._table.item(row, 1).text()
            src_combo = self._table.cellWidget(row, 2)
            default_combo = self._table.cellWidget(row, 3)
            if src_combo and default_combo and src_combo.currentIndex() == 0:
                result[ptype] = default_combo.currentText()
        return result

    def get_selected_files(self) -> List[str]:
        return [f for f, cb in self._file_checks.items() if cb.isChecked()]


# ══════════════════════════════════════════════════════════════
#   File Editor Tab
# ══════════════════════════════════════════════════════════════

class FileEditorTab(QWidget):
    """A single tab in the file editor showing one OpenFOAM file."""

    saved = pyqtSignal(str)   # emits file path when saved
    modified_changed = pyqtSignal(bool)

    def __init__(self, file_path: str, parent=None):
        super().__init__(parent)
        self.file_path = file_path
        self._modified = False
        self._build_ui()
        self._load()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        # Toolbar
        bar = QToolBar()
        bar.setIconSize(QSize(16, 16))
        bar.setStyleSheet('QToolBar { border:none; background:#2d2d30; padding:2px; }')

        self._path_label = QLabel(f'  {self.file_path}')
        self._path_label.setStyleSheet('color:#9cdcfe; font-size:11px;')
        bar.addWidget(self._path_label)

        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        bar.addWidget(spacer)

        save_btn = QPushButton('Save')
        save_btn.setFixedHeight(24)
        save_btn.setStyleSheet(
            'QPushButton { background:#0e639c; color:white; border:none; '
            'padding:2px 10px; border-radius:3px; }'
            'QPushButton:hover { background:#1177bb; }'
        )
        save_btn.clicked.connect(self._save)
        bar.addWidget(save_btn)

        reload_btn = QPushButton('Reload')
        reload_btn.setFixedHeight(24)
        reload_btn.setStyleSheet(
            'QPushButton { background:#3c3c3c; color:#d4d4d4; border:none; '
            'padding:2px 10px; border-radius:3px; }'
            'QPushButton:hover { background:#555; }'
        )
        reload_btn.clicked.connect(self._reload)
        bar.addWidget(reload_btn)

        layout.addWidget(bar)

        # Editor
        self._editor = CodeEditor()
        self._editor.document().contentsChanged.connect(self._mark_modified)
        layout.addWidget(self._editor)

    def _load(self):
        try:
            text = Path(self.file_path).read_text(encoding='utf-8', errors='replace')
            self._editor.document().contentsChanged.disconnect(self._mark_modified)
            self._editor.setPlainText(text)
            self._editor.document().contentsChanged.connect(self._mark_modified)
            self._modified = False
            self.modified_changed.emit(False)
        except Exception as exc:
            self._editor.setPlainText(f'// Error loading file:\n// {exc}')

    def _save(self):
        try:
            Path(self.file_path).write_text(
                self._editor.toPlainText(), encoding='utf-8'
            )
            self._modified = False
            self.modified_changed.emit(False)
            self.saved.emit(self.file_path)
        except Exception as exc:
            QMessageBox.critical(self, 'Save Error', str(exc))

    def _reload(self):
        if self._modified:
            ans = QMessageBox.question(
                self, 'Discard Changes',
                'Reload file and discard unsaved changes?',
                QMessageBox.Yes | QMessageBox.No,
            )
            if ans != QMessageBox.Yes:
                return
        self._load()

    def _mark_modified(self):
        if not self._modified:
            self._modified = True
            self.modified_changed.emit(True)

    def reload_content(self):
        """Public method to reload the file (e.g. after BC update)."""
        self._load()

    @property
    def is_modified(self) -> bool:
        return self._modified


# ══════════════════════════════════════════════════════════════
#   Editor Panel (tab widget holding multiple file editors)
# ══════════════════════════════════════════════════════════════

class EditorPanel(QTabWidget):
    status_message = pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setTabsClosable(True)
        self.setMovable(True)
        self.setDocumentMode(True)
        self.tabCloseRequested.connect(self._close_tab)
        self._open_files: Dict[str, int] = {}  # path -> tab index

        self.setStyleSheet(
            'QTabWidget::pane { border-top: 1px solid #3e3e42; }'
            'QTabBar::tab { background:#2d2d30; color:#9d9d9d; '
            '  padding:5px 14px; border-right:1px solid #3e3e42; }'
            'QTabBar::tab:selected { background:#1e1e1e; color:#fff; '
            '  border-top:1px solid #007acc; }'
            'QTabBar::tab:hover { background:#383838; }'
        )

    def open_file(self, file_path: str):
        """Open a file in the editor (or focus existing tab)."""
        abs_path = str(Path(file_path).resolve())
        if abs_path in self._open_files:
            self.setCurrentIndex(self._open_files[abs_path])
            return

        tab = FileEditorTab(abs_path)
        tab.saved.connect(lambda p: self.status_message.emit(f'Saved: {p}'))
        tab.modified_changed.connect(
            lambda mod, p=abs_path: self._update_tab_title(p, mod)
        )
        idx = self.addTab(tab, Path(abs_path).name)
        self._open_files[abs_path] = idx
        self.setCurrentIndex(idx)
        self.setTabToolTip(idx, abs_path)

    def reload_file(self, file_path: str):
        """Reload a file if it is currently open."""
        abs_path = str(Path(file_path).resolve())
        if abs_path in self._open_files:
            tab = self.widget(self._open_files[abs_path])
            if isinstance(tab, FileEditorTab):
                tab.reload_content()

    def _close_tab(self, idx: int):
        tab = self.widget(idx)
        if isinstance(tab, FileEditorTab) and tab.is_modified:
            ans = QMessageBox.question(
                self, 'Unsaved Changes',
                f'Close "{Path(tab.file_path).name}" without saving?',
                QMessageBox.Yes | QMessageBox.No,
            )
            if ans != QMessageBox.Yes:
                return
        # Remove from tracking dict
        path_to_remove = None
        for p, i in self._open_files.items():
            if i == idx:
                path_to_remove = p
                break
        if path_to_remove:
            del self._open_files[path_to_remove]
        # Update indices for tabs after the removed one
        self._open_files = {
            p: (i - 1 if i > idx else i)
            for p, i in self._open_files.items()
        }
        self.removeTab(idx)

    def _update_tab_title(self, path: str, modified: bool):
        if path not in self._open_files:
            return
        idx = self._open_files[path]
        base = Path(path).name
        self.setTabText(idx, ('● ' if modified else '') + base)


# ══════════════════════════════════════════════════════════════
#   Case Structure Tree
# ══════════════════════════════════════════════════════════════

class CaseTreeWidget(QTreeWidget):
    file_activated = pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setHeaderLabel('Case Files')
        self.setAnimated(True)
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self._context_menu)
        self.itemDoubleClicked.connect(self._on_double_click)
        self.setStyleSheet(
            'QTreeWidget { background:#252526; color:#ccc; border:none; }'
            'QTreeWidget::item:hover { background:#2a2d2e; }'
            'QTreeWidget::item:selected { background:#094771; }'
        )
        self._case_path = ''

    def load_case(self, case_path: str):
        self._case_path = case_path
        self.clear()
        p = Path(case_path)
        root = QTreeWidgetItem([p.name])
        root.setData(0, Qt.UserRole, str(p))
        self.addTopLevelItem(root)
        self._populate(root, p)
        root.setExpanded(True)

    def _populate(self, parent_item: QTreeWidgetItem, directory: Path, depth=0):
        if depth > 4:
            return
        try:
            for item in sorted(directory.iterdir(), key=lambda x: (x.is_file(), x.name)):
                child = QTreeWidgetItem([item.name])
                child.setData(0, Qt.UserRole, str(item))
                if item.is_dir():
                    child.setForeground(0, QColor('#c5a5c5'))
                    self._populate(child, item, depth + 1)
                else:
                    child.setForeground(0, QColor('#9cdcfe'))
                parent_item.addChild(child)
        except PermissionError:
            pass

    def _on_double_click(self, item: QTreeWidgetItem, col: int):
        path = item.data(0, Qt.UserRole)
        if path and Path(path).is_file():
            self.file_activated.emit(path)

    def _context_menu(self, pos: QPoint):
        item = self.itemAt(pos)
        if not item:
            return
        path = item.data(0, Qt.UserRole)
        if not path or not Path(path).is_file():
            return
        menu = QMenu(self)
        open_act = menu.addAction('Open in Editor')
        act = menu.exec_(self.viewport().mapToGlobal(pos))
        if act == open_act:
            self.file_activated.emit(path)


# ══════════════════════════════════════════════════════════════
#   Left Control Panel
# ══════════════════════════════════════════════════════════════

class ControlPanel(QWidget):
    open_file_requested = pyqtSignal(str)
    status_message = pyqtSignal(str)
    case_changed = pyqtSignal(str, str)  # source, target

    def __init__(self, editor: EditorPanel, parent=None):
        super().__init__(parent)
        self._editor = editor
        self._source = ''
        self._target = ''
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(8)

        # ── Source case ───────────────────────────────────────
        src_group = QGroupBox('Source Case')
        src_layout = QVBoxLayout(src_group)

        src_row = QHBoxLayout()
        self._src_edit = QLineEdit()
        self._src_edit.setPlaceholderText('Path to existing OpenFOAM case...')
        self._src_edit.editingFinished.connect(self._on_src_changed)
        src_browse = QPushButton('Browse')
        src_browse.setFixedWidth(70)
        src_browse.clicked.connect(self._browse_source)
        src_row.addWidget(self._src_edit)
        src_row.addWidget(src_browse)
        src_layout.addLayout(src_row)

        self._src_info = QLabel('<i>No case loaded</i>')
        self._src_info.setWordWrap(True)
        self._src_info.setStyleSheet('color:#9cdcfe; font-size:11px;')
        src_layout.addWidget(self._src_info)
        layout.addWidget(src_group)

        # ── Target case ───────────────────────────────────────
        tgt_group = QGroupBox('Target Case (New Project)')
        tgt_layout = QVBoxLayout(tgt_group)

        tgt_row = QHBoxLayout()
        self._tgt_edit = QLineEdit()
        self._tgt_edit.setPlaceholderText('Path to new/target OpenFOAM case...')
        self._tgt_edit.editingFinished.connect(self._on_tgt_changed)
        tgt_browse = QPushButton('Browse')
        tgt_browse.setFixedWidth(70)
        tgt_browse.clicked.connect(self._browse_target)
        tgt_row.addWidget(self._tgt_edit)
        tgt_row.addWidget(tgt_browse)
        tgt_layout.addLayout(tgt_row)

        self._tgt_info = QLabel('<i>No case loaded</i>')
        self._tgt_info.setWordWrap(True)
        self._tgt_info.setStyleSheet('color:#9cdcfe; font-size:11px;')
        tgt_layout.addWidget(self._tgt_info)
        layout.addWidget(tgt_group)

        # ── Action buttons ────────────────────────────────────
        act_group = QGroupBox('Actions')
        act_layout = QVBoxLayout(act_group)

        self._copy_btn = QPushButton('1.  Copy Case Settings  →  Target')
        self._copy_btn.setToolTip(
            'Copy system/, 0/, and constant/ (excluding polyMesh) from source to target'
        )
        self._copy_btn.setStyleSheet(
            'QPushButton { background:#0e639c; color:white; padding:6px; '
            'border-radius:3px; font-weight:bold; }'
            'QPushButton:hover { background:#1177bb; }'
            'QPushButton:disabled { background:#3c3c3c; color:#666; }'
        )
        self._copy_btn.setEnabled(False)
        self._copy_btn.clicked.connect(self._do_copy)
        act_layout.addWidget(self._copy_btn)

        self._map_btn = QPushButton('2.  Map Boundary Conditions')
        self._map_btn.setToolTip(
            'Open dialog to map source patches to target mesh patches and update BC files'
        )
        self._map_btn.setStyleSheet(
            'QPushButton { background:#5c7a2e; color:white; padding:6px; '
            'border-radius:3px; font-weight:bold; }'
            'QPushButton:hover { background:#6d9236; }'
            'QPushButton:disabled { background:#3c3c3c; color:#666; }'
        )
        self._map_btn.setEnabled(False)
        self._map_btn.clicked.connect(self._do_map)
        act_layout.addWidget(self._map_btn)

        open_tgt_btn = QPushButton('Open Target File...')
        open_tgt_btn.setStyleSheet(
            'QPushButton { background:#3c3c3c; color:#d4d4d4; padding:5px; '
            'border-radius:3px; }'
            'QPushButton:hover { background:#555; }'
        )
        open_tgt_btn.clicked.connect(self._open_file)
        act_layout.addWidget(open_tgt_btn)

        layout.addWidget(act_group)

        # ── Case tree ─────────────────────────────────────────
        tree_label = QLabel('Target Case Structure:')
        tree_label.setStyleSheet('color:#858585; font-size:11px; margin-top:4px;')
        layout.addWidget(tree_label)

        self._tree = CaseTreeWidget()
        self._tree.file_activated.connect(self._editor.open_file)
        layout.addWidget(self._tree)

    # ── Slots ──────────────────────────────────────────────────

    def _browse_source(self):
        d = QFileDialog.getExistingDirectory(
            self, 'Select Source OpenFOAM Case',
            self._source or str(Path.home()),
        )
        if d:
            self._src_edit.setText(d)
            self._on_src_changed()

    def _browse_target(self):
        d = QFileDialog.getExistingDirectory(
            self, 'Select Target OpenFOAM Case Directory',
            self._target or self._source or str(Path.home()),
        )
        if d:
            self._tgt_edit.setText(d)
            self._on_tgt_changed()

    def _on_src_changed(self):
        path = self._src_edit.text().strip()
        if not path:
            return
        if not Path(path).is_dir():
            self._src_info.setText('<span style="color:#f44">Directory not found.</span>')
            return
        self._source = path
        info = OpenFOAMParser.get_case_info(path)
        if not info and not OpenFOAMParser.is_valid_case(path):
            self._src_info.setText(
                '<span style="color:#fa0">⚠ Does not look like an OpenFOAM case '
                '(missing system/ or constant/).</span>'
            )
        else:
            parts = []
            if 'application' in info:
                parts.append(f"app: <b>{info['application']}</b>")
            if 'patch_count' in info:
                parts.append(f"patches: <b>{info['patch_count']}</b> ({info.get('patches', '')})")
            self._src_info.setText('  '.join(parts) or '<i>Valid case</i>')
        self._update_buttons()
        self.case_changed.emit(self._source, self._target)

    def _on_tgt_changed(self):
        path = self._tgt_edit.text().strip()
        if not path:
            return
        self._target = path
        p = Path(path)
        if p.is_dir():
            info = OpenFOAMParser.get_case_info(path)
            if info:
                patches_text = info.get('patches', '')
                self._tgt_info.setText(
                    f"patches: <b>{info.get('patch_count', '?')}</b> ({patches_text})"
                )
            else:
                self._tgt_info.setText('<i>Directory exists (no case yet)</i>')
            self._tree.load_case(path)
        else:
            self._tgt_info.setText('<i>Will be created on copy</i>')
        self._update_buttons()
        self.case_changed.emit(self._source, self._target)

    def _update_buttons(self):
        has_src = bool(self._source) and Path(self._source).is_dir()
        has_tgt = bool(self._target)
        self._copy_btn.setEnabled(has_src and has_tgt)

        has_tgt_mesh = (
            has_tgt
            and (Path(self._target) / 'constant' / 'polyMesh' / 'boundary').exists()
        )
        self._map_btn.setEnabled(has_src and has_tgt and has_tgt_mesh)

    def _do_copy(self):
        if not self._source or not self._target:
            return
        if self._source == self._target:
            QMessageBox.warning(self, 'Same Path', 'Source and target must be different.')
            return

        ans = QMessageBox.question(
            self, 'Confirm Copy',
            f'Copy case settings from:\n  {self._source}\nto:\n  {self._target}\n\n'
            f'(existing system/, 0/, 0.orig/, and non-polyMesh constant/ files '
            f'in the target will be overwritten)',
            QMessageBox.Yes | QMessageBox.No,
        )
        if ans != QMessageBox.Yes:
            return

        self._copy_btn.setEnabled(False)
        self._progress = QProgressDialog('Copying case files...', None, 0, 0, self)
        self._progress.setWindowTitle('Copying...')
        self._progress.setWindowModality(Qt.WindowModal)
        self._progress.show()

        self._worker = CopyWorker(self._source, self._target)
        self._worker.progress.connect(
            lambda msg: self._progress.setLabelText(msg)
        )
        self._worker.finished.connect(self._on_copy_done)
        self._worker.start()

    def _on_copy_done(self, success: bool, message: str):
        self._progress.close()
        self._copy_btn.setEnabled(True)
        if success:
            self.status_message.emit(message)
            self._tree.load_case(self._target)
            self._update_buttons()
            QMessageBox.information(self, 'Copy Complete', message)
        else:
            QMessageBox.critical(self, 'Copy Failed', message)

    def _do_map(self):
        if not self._source or not self._target:
            return

        src_patches = OpenFOAMParser.get_boundary_patches(self._source)
        tgt_patches = OpenFOAMParser.get_boundary_patches(self._target)

        if not tgt_patches:
            QMessageBox.warning(
                self, 'No Target Mesh',
                'Could not read patches from target constant/polyMesh/boundary.\n'
                'Make sure the new mesh is placed in the target case first.',
            )
            return

        bc_files, zero_dir = OpenFOAMParser.get_bc_files(self._target)
        if not bc_files:
            QMessageBox.warning(
                self, 'No BC Files',
                'No field files found in target 0/ directory.\n'
                'Copy the case settings first.',
            )
            return

        dlg = BoundaryMappingDialog(
            self._source, self._target,
            src_patches, tgt_patches, bc_files, self
        )
        if dlg.exec_() != QDialog.Accepted:
            return

        mapping = dlg.get_mapping()
        default_bc_by_type = dlg.get_default_bc_by_type()
        selected_files = dlg.get_selected_files()

        if not selected_files:
            QMessageBox.information(self, 'Nothing to do', 'No files were selected.')
            return

        errors = []
        for fname in selected_files:
            fpath = os.path.join(zero_dir, fname)
            ok, msg = OpenFOAMParser.update_boundary_field(
                fpath, mapping, tgt_patches, default_bc_by_type
            )
            if not ok:
                errors.append(f'{fname}: {msg}')
            else:
                # Refresh in editor if open
                self._editor.reload_file(fpath)

        self._tree.load_case(self._target)

        if errors:
            QMessageBox.warning(
                self, 'Partial Errors',
                'Some files could not be updated:\n' + '\n'.join(errors),
            )
            self.status_message.emit(f'BC mapping done with {len(errors)} error(s).')
        else:
            self.status_message.emit(
                f'Boundary conditions updated in {len(selected_files)} file(s).'
            )
            QMessageBox.information(
                self, 'Done',
                f'Boundary conditions updated in {len(selected_files)} file(s).',
            )

    def _open_file(self):
        start_dir = self._target if self._target else str(Path.home())
        path, _ = QFileDialog.getOpenFileName(
            self, 'Open OpenFOAM File', start_dir,
            'All Files (*)',
        )
        if path:
            self._editor.open_file(path)

    def refresh_tree(self):
        if self._target and Path(self._target).is_dir():
            self._tree.load_case(self._target)


# ══════════════════════════════════════════════════════════════
#   Main Window
# ══════════════════════════════════════════════════════════════

DARK_STYLE = """
QMainWindow, QDialog, QWidget {
    background: #1e1e1e;
    color: #d4d4d4;
}
QGroupBox {
    border: 1px solid #3e3e42;
    border-radius: 4px;
    margin-top: 10px;
    padding-top: 4px;
    font-weight: bold;
    color: #d4d4d4;
}
QGroupBox::title {
    subcontrol-origin: margin;
    left: 8px;
    padding: 0 4px;
    color: #9cdcfe;
}
QLineEdit, QComboBox, QPlainTextEdit {
    background: #3c3c3c;
    color: #d4d4d4;
    border: 1px solid #555;
    border-radius: 3px;
    padding: 3px;
    selection-background-color: #094771;
}
QLineEdit:focus, QComboBox:focus {
    border: 1px solid #007acc;
}
QPushButton {
    background: #3c3c3c;
    color: #d4d4d4;
    border: 1px solid #555;
    border-radius: 3px;
    padding: 4px 8px;
}
QPushButton:hover { background: #555; }
QPushButton:pressed { background: #094771; }
QScrollBar:vertical {
    background: #252526;
    width: 12px;
}
QScrollBar::handle:vertical {
    background: #424242;
    border-radius: 4px;
    min-height: 20px;
}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0; }
QScrollBar:horizontal {
    background: #252526;
    height: 12px;
}
QScrollBar::handle:horizontal {
    background: #424242;
    border-radius: 4px;
    min-width: 20px;
}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal { width: 0; }
QTreeWidget { background: #252526; }
QHeaderView::section {
    background: #252526;
    color: #9cdcfe;
    padding: 4px;
    border: none;
    border-bottom: 1px solid #3e3e42;
}
QTableWidget {
    background: #252526;
    gridline-color: #3e3e42;
    selection-background-color: #094771;
}
QTableWidget::item:alternate { background: #2a2d2e; }
QStatusBar { background: #007acc; color: white; }
QStatusBar QLabel { color: white; padding: 0 6px; }
QSplitter::handle { background: #3e3e42; }
QCheckBox { spacing: 6px; }
QCheckBox::indicator {
    width: 14px; height: 14px;
    border: 1px solid #555;
    border-radius: 2px;
    background: #3c3c3c;
}
QCheckBox::indicator:checked { background: #007acc; border-color: #007acc; }
QMenu {
    background: #252526;
    border: 1px solid #454545;
}
QMenu::item:selected { background: #094771; }
QProgressDialog { background: #252526; }
QDialogButtonBox QPushButton { min-width: 80px; }
"""


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('OpenFOAM Case Manager')
        self.setMinimumSize(1200, 720)
        self.resize(1440, 820)
        self._build_ui()
        self._build_menu()
        self._status_bar = self.statusBar()
        self._status_bar.showMessage('Ready — select a source and target case to begin.')

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        splitter = QSplitter(Qt.Horizontal)
        splitter.setHandleWidth(4)

        # Editor panel (right side, created first so ControlPanel can ref it)
        self._editor = EditorPanel()
        self._editor.status_message.connect(
            lambda m: self.statusBar().showMessage(m, 5000)
        )

        # Control panel (left side)
        self._ctrl = ControlPanel(self._editor)
        self._ctrl.status_message.connect(
            lambda m: self.statusBar().showMessage(m, 6000)
        )
        self._ctrl.setMinimumWidth(300)
        self._ctrl.setMaximumWidth(500)

        # Welcome / placeholder for right side
        placeholder = QWidget()
        ph_layout = QVBoxLayout(placeholder)
        ph_layout.setAlignment(Qt.AlignCenter)
        welcome = QLabel(
            '<h2 style="color:#569cd6;">OpenFOAM Case Manager</h2>'
            '<p style="color:#9cdcfe;">Quick-copy simulation settings between cases<br>'
            'and adapt boundary conditions for new meshes.</p>'
            '<p style="color:#858585;">1. Set <b>Source Case</b> (existing case)<br>'
            '2. Set <b>Target Case</b> (new case directory)<br>'
            '3. Click <b>Copy Case Settings</b><br>'
            '4. Click <b>Map Boundary Conditions</b><br>'
            '5. Edit files directly in this editor</p>'
        )
        welcome.setAlignment(Qt.AlignCenter)
        ph_layout.addWidget(welcome)

        # Tab widget with welcome and editor
        self._tabs = QTabWidget()
        self._tabs.setTabPosition(QTabWidget.North)
        self._tabs.setDocumentMode(True)
        self._tabs.addTab(placeholder, 'Welcome')
        self._tabs.addTab(self._editor, 'File Editor')
        self._tabs.setStyleSheet(
            'QTabWidget::pane { border-top: 1px solid #3e3e42; }'
            'QTabBar::tab { background:#2d2d30; color:#9d9d9d; '
            '  padding:6px 16px; }'
            'QTabBar::tab:selected { background:#1e1e1e; color:#fff; '
            '  border-top:1px solid #007acc; }'
        )

        # When a file is opened in the editor, switch to editor tab
        self._editor.currentChanged.connect(self._on_editor_tab_change)

        splitter.addWidget(self._ctrl)
        splitter.addWidget(self._tabs)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([340, 1060])

        root.addWidget(splitter)

    def _on_editor_tab_change(self, idx: int):
        if idx >= 0:
            self._tabs.setCurrentIndex(1)  # switch to File Editor tab

    def _build_menu(self):
        mb = self.menuBar()
        mb.setStyleSheet(
            'QMenuBar { background:#2d2d30; color:#d4d4d4; }'
            'QMenuBar::item:selected { background:#094771; }'
            'QMenu { background:#252526; border:1px solid #454545; }'
            'QMenu::item:selected { background:#094771; }'
        )

        file_menu = mb.addMenu('&File')

        open_act = QAction('Open File...', self)
        open_act.setShortcut('Ctrl+O')
        open_act.triggered.connect(self._menu_open_file)
        file_menu.addAction(open_act)

        file_menu.addSeparator()
        quit_act = QAction('Quit', self)
        quit_act.setShortcut('Ctrl+Q')
        quit_act.triggered.connect(self.close)
        file_menu.addAction(quit_act)

        help_menu = mb.addMenu('&Help')
        about_act = QAction('About', self)
        about_act.triggered.connect(self._show_about)
        help_menu.addAction(about_act)

    def _menu_open_file(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Open File', '', 'All Files (*)')
        if path:
            self._editor.open_file(path)
            self._tabs.setCurrentIndex(1)

    def _show_about(self):
        QMessageBox.about(
            self,
            'About OpenFOAM Case Manager',
            '<h3>OpenFOAM Case Manager</h3>'
            '<p>A PyQt5-based tool for managing OpenFOAM simulation cases.</p>'
            '<ul>'
            '<li>Copy case settings (excluding polyMesh)</li>'
            '<li>Map boundary conditions to new mesh patches</li>'
            '<li>Edit OpenFOAM dictionaries with syntax highlighting</li>'
            '</ul>'
            '<p style="color:#858;">Supports OpenFOAM v5+ / ESI-OpenFOAM</p>',
        )


# ══════════════════════════════════════════════════════════════
#   Entry Point
# ══════════════════════════════════════════════════════════════

def main():
    app = QApplication(sys.argv)
    app.setApplicationName('OpenFOAM Case Manager')
    app.setStyle(QStyleFactory.create('Fusion'))

    # Apply dark palette
    palette = QPalette()
    palette.setColor(QPalette.Window, QColor('#1e1e1e'))
    palette.setColor(QPalette.WindowText, QColor('#d4d4d4'))
    palette.setColor(QPalette.Base, QColor('#252526'))
    palette.setColor(QPalette.AlternateBase, QColor('#2d2d30'))
    palette.setColor(QPalette.ToolTipBase, QColor('#252526'))
    palette.setColor(QPalette.ToolTipText, QColor('#d4d4d4'))
    palette.setColor(QPalette.Text, QColor('#d4d4d4'))
    palette.setColor(QPalette.Button, QColor('#3c3c3c'))
    palette.setColor(QPalette.ButtonText, QColor('#d4d4d4'))
    palette.setColor(QPalette.BrightText, QColor('#ffffff'))
    palette.setColor(QPalette.Link, QColor('#007acc'))
    palette.setColor(QPalette.Highlight, QColor('#094771'))
    palette.setColor(QPalette.HighlightedText, QColor('#ffffff'))
    app.setPalette(palette)
    app.setStyleSheet(DARK_STYLE)

    window = MainWindow()
    window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
