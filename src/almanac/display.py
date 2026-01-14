import logging
import numpy as np
from itertools import cycle
from datetime import datetime, timedelta
from rich.console import Console
from rich.table import Table
from rich.live import Live
from rich.text import Text
from rich.align import Align
from rich.table import Table as RichTable

from almanac import config
from almanac.data_models import Exposure
from typing import Optional, List, Tuple, Dict, Any

def mjd_to_datetime(mjd):
    """Convert MJD to datetime - mock implementation"""
    # MJD epoch is November 17, 1858
    epoch = datetime(1858, 11, 17)
    return epoch + timedelta(days=mjd)

class BufferedHandler(logging.Handler):
    """Custom logging handler that buffers log records"""
    def __init__(self):
        super().__init__()
        self.buffer = []

    def emit(self, record):
        self.buffer.append(record)

    def flush_to_console(self, console=None):
        """Flush buffered records to console"""
        if console is None:
            console = Console()

        for record in self.buffer:
            log_message = self.format(record)

            # Color code based on log level
            if record.levelno >= logging.ERROR:
                style = "red"
            elif record.levelno >= logging.WARNING:
                style = "yellow"
            elif record.levelno >= logging.INFO:
                style = "blue"
            else:
                style = "dim"

            console.print(log_message, style=style)

        self.buffer.clear()

class ObservationsDisplay:
    color_outside_range = "black"
    color_unknown = "white"
    color_no_data = "bright_black"
    color_apo = "dodger_blue3"
    color_lco = "green4"
    color_both = "purple4"
    color_missing = "red"

    def __init__(self, mjd_min, mjd_max, observatories=("apo", "lco")):
        self.console = Console()
        self.start_date = mjd_to_datetime(mjd_min) if isinstance(mjd_min, (int, float)) else mjd_min
        self.end_date = mjd_to_datetime(mjd_max) if isinstance(mjd_max, (int, float)) else mjd_max
        self.days_per_week = 7

        # Track completion status for each day
        self.completed = dict(apo=set(), lco=set())
        self.no_data = dict(apo=set(), lco=set())
        self.missing = set()
        self.offset = 0
        # Setup logging buffer
        self.log_buffer = BufferedHandler()
        self.log_buffer.setFormatter(logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%H:%M:%S'
        ))

        # Calculate the grid dimensions based on the date range
        self._setup_grid()
        self.observatories = observatories

    def _setup_grid(self):
        """Setup the grid based on start and end dates, organized by years"""
        # Calculate which years we need to display
        self.start_year = self.start_date.year
        self.end_year = self.end_date.year
        self.year_rows = list(range(self.start_year, self.end_year + 1))

        # For each year, we'll create a full year grid (Jan 1 to Dec 31)
        # but pad the beginning if our data doesn't start on Jan 1
        self.year_grids = {}
        self.dates = []

        self.offset = 0
        for i, year in enumerate(self.year_rows):
            year_start = datetime(year, 1, 1)
            year_end = datetime(year, 12, 31)

            # Find Monday of the week containing Jan 1
            jan1_weekday = year_start.weekday()  # 0=Monday, 6=Sunday
            grid_start = year_start - timedelta(days=jan1_weekday)
            if i == 0:
                self.offset = (self.start_date - grid_start).days

            # Find Sunday of the week containing Dec 31
            dec31_weekday = year_end.weekday()  # 0=Monday, 6=Sunday
            days_to_sunday = 6 - dec31_weekday
            grid_end = year_end + timedelta(days=days_to_sunday)

            # Generate dates for this year's grid
            year_dates = []
            current_date = grid_start
            while current_date <= grid_end:
                year_dates.append(current_date)
                current_date += timedelta(days=1)

            # Calculate weeks for this year
            total_days = len(year_dates)
            weeks = total_days // 7

            self.year_grids[year] = {
                'dates': year_dates,
                'weeks': weeks,
                'grid_start': grid_start,
                'grid_end': grid_end,
                'year_start': year_start,
                'year_end': year_end
            }

            # Add to master dates list with offset
            for i, date in enumerate(year_dates):
                self.dates.append(date)


    def get_day_color(self, day_index):
        """Return the color for a given day based on completion status"""
        if day_index >= len(self.dates):
            return self.color_outside_range

        date = self.dates[day_index]

        # Only show colored squares for dates within our actual range
        if date < self.start_date or date > self.end_date:
            return self.color_outside_range

        if day_index in self.missing:
            return self.color_missing
        elif day_index in self.completed["apo"] and day_index in self.completed["lco"]:
            return self.color_both
        elif day_index in self.completed["apo"]:
            return self.color_apo
        elif day_index in self.completed["lco"]:
            return self.color_lco
        elif day_index in self.no_data["apo"] and day_index in self.no_data["lco"]:
            return self.color_no_data
        else:
            return self.color_unknown

    def _get_month_headers_for_year(self, year):
        """Generate month headers for a specific year"""
        year_data = self.year_grids[year]
        weeks = year_data['weeks']
        year_dates = year_data['dates']

        headers = ["   "]  # Space for day labels
        current_month = None
        text_to_add = ""

        for week in range(weeks):
            week_start_index = week * 7

            if week_start_index < len(year_dates):
                week_date = year_dates[week_start_index]

                # Only show month headers for dates within the actual year
                if week_date.year == year:
                    month_abbr = week_date.strftime("%b")

                    # Only show month if it's different from previous week
                    if current_month != month_abbr:
                        headers.append(Text(f"{month_abbr[:1]}", style="dim"))
                        text_to_add = month_abbr[1:]
                        current_month = month_abbr
                    else:
                        if len(text_to_add) > 0:
                            headers.append(Text(text_to_add[:1], style="dim"))
                            text_to_add = text_to_add[1:]
                        else:
                            headers.append(Text(" "))
                else:
                    # This is padding (before Jan 1 or after Dec 31)
                    headers.append(Text(" "))
            else:
                headers.append(Text(" "))

        return headers

    def create_contributions_grid_for_year(self, year):
        """Create the contributions grid for a specific year"""
        year_data = self.year_grids[year]
        weeks = year_data['weeks']
        year_dates = year_data['dates']

        table = Table.grid(padding=0)

        # Add columns for day labels and each week in this year
        table.add_column()  # For day labels
        for _ in range(weeks):
            table.add_column()

        # Add month headers
        month_headers = self._get_month_headers_for_year(year)
        table.add_row(*month_headers)

        # Create rows for each day of the week
        day_names = ["S", "M", "T", "W", "T", "F", "S"]

        for day_of_week in range(self.days_per_week):
            row = [Text(day_names[day_of_week].ljust(3), style="dim")]

            for week in range(weeks):
                day_index_in_year = week * 7 + day_of_week

                if day_index_in_year < len(year_dates):
                    date = year_dates[day_index_in_year]

                    # Find this date in our master dates list to get the right index
                    master_day_index = None
                    for i, master_date in enumerate(self.dates):
                        if master_date == date:
                            master_day_index = i
                            break

                    if master_day_index is not None:
                        color = self.get_day_color(master_day_index)

                        # Show square only for dates within the actual year and our date range
                        if (date.year == year and
                            self.start_date <= date <= self.end_date):
                            square = Text("■", style=color)
                        elif date.year == year:
                            # Within the year but outside our date range
                            square = Text("■", style=self.color_no_data)
                        else:
                            # Padding dates (before Jan 1 or after Dec 31)
                            square = Text(" ")
                    else:
                        square = Text(" ")
                else:
                    square = Text(" ")

                row.append(square)

            table.add_row(*row)

        return table

    def create_display(self):
        """Create the complete display with title and yearly grids"""
        date_range = f"{self.start_date.strftime('%b %d, %Y')} - {self.end_date.strftime('%b %d, %Y')}"
        title = Text("SDSS/APOGEE Observations", style="bold white")
        subtitle = Text(date_range, style="dim")

        # Create legend
        legend = Table.grid(padding=(0, 1))
        legend.add_column()
        legend.add_column()
        legend.add_column()
        legend.add_column()
        legend.add_column()

        items = [
            Text("■", style=self.color_no_data),
            Text("None", style="dim"),
        ]
        if "apo" in self.observatories:
            items.extend([
                Text("■", style=self.color_apo),
                Text("APO", style="dim"),
            ])
        if "lco" in self.observatories:
            items.extend([
                Text("■", style=self.color_lco),
                Text("LCO", style="dim"),
            ])

        if "apo" in self.observatories and "lco" in self.observatories:
            items.extend([
                Text("■", style=self.color_both),
                Text("Both", style="dim"),
            ])
        legend.add_row(*items)

        # Combine everything
        main_table = Table.grid()
        main_table.add_column()
        main_table.add_row(Align.center(title))
        main_table.add_row(Align.center(subtitle))
        main_table.add_row("")

        # Add each year's grid with year header
        for i, year in enumerate(self.year_rows):
            # Add year header
            year_header = Text(str(year), style="bold cyan")
            main_table.add_row(Align.left(year_header))

            # Add the grid for this year
            year_grid = self.create_contributions_grid_for_year(year)
            main_table.add_row(Align.center(year_grid))

            # Add spacing between years except after the last one
            if i < len(self.year_rows) - 1:
                main_table.add_row("")

        main_table.add_row("")
        main_table.add_row(Align.center(legend))

        return main_table

    def add_observation(self, date, observatory):
        """Add an observation for a specific date and observatory"""
        # Find the day index for this date
        for i, grid_date in enumerate(self.dates):
            if grid_date.date() == date.date():
                self.completed[observatory].add(i)
                break


def display_exposures(
    exposures: List[Exposure],
    sequences: Optional[Dict[str, List[Tuple[int, int]]]] = None,
    console: Optional[Console] = None,
    header_style: str = "bold cyan",
    column_names: Optional[List[str]] = None,
    sequence_styles: Tuple[str, ...] = ("green", "yellow"),
    missing_style: str = "blink bold red",
    title_style: str = "bold blue",
) -> None:
    """Display exposure information using Rich table formatting.

    Args:
        exposures: List of Exposure objects containing exposure data
        sequences: Dictionary mapping sequence names to lists of (start, end) tuples (default: None)
        console: Rich Console instance (default: None, creates new one)
        header_style: Style for table headers (default: "bold cyan")
        sequence_styles: Tuple of styles to cycle through for sequences (default: ("green", "yellow"))
        missing_style: Style for missing/error entries (default: "red")
        title_style: Style for the table title (default: "bold blue")
    """
    if console is None:
        console = Console()

    if len(exposures) == 0:
        return

    # Create the title
    observatory, mjd = (exposures[0].observatory, exposures[0].mjd)
    title = f"{len(exposures)} exposures from {observatory.upper()} on MJD {mjd}"

    # Create Rich table
    rich_table = RichTable(title=title, title_style=title_style, show_header=True, header_style=header_style)

    field_names = config.display_field_names

    for field_name in field_names:
        rich_table.add_column(field_name, justify="center")

    # Prepare sequence tracking
    flattened_sequences = []
    for k, v in (sequences or dict()).items():
        flattened_sequences.extend(v)
    flattened_sequences = np.array(flattened_sequences)

    sequence_styles_cycle = cycle(sequence_styles)
    in_sequence, current_sequence_style = (False, next(sequence_styles_cycle))

    # Add rows to the table
    for i, exposure in enumerate(exposures, start=1):
        # Check if this row is part of a sequence
        row_style = None
        end_of_sequence = None
        if len(flattened_sequences) > 0:
            try:
                j, k = np.where(flattened_sequences == i)
            except:
                pass
            else:
                # Could be start or end of sequence, and could be out of order
                start_of_sequence = 0 in k
                end_of_sequence = 1 in k

                if start_of_sequence:
                    in_sequence = True
                    current_sequence_style = next(sequence_styles_cycle)
                elif end_of_sequence:  # only end of sequence
                    in_sequence = True

        # Determine row style
        if in_sequence:
            row_style = current_sequence_style
        else:
            # Check if it's missing or has issues
            if exposure.image_type == "missing":
                row_style = missing_style

        # Convert row data to strings and apply styling if needed
        row_data = []
        for field_name in field_names:
            value = getattr(exposure, field_name)
            if row_style:
                row_data.append(Text(f"{value}", style=row_style))
            else:
                row_data.append(f"{value}")

        rich_table.add_row(*row_data)
        if end_of_sequence:
            in_sequence = False

    console.print(rich_table)
    console.print()  # Add a blank line after the table


from rich.live import Live
from rich.table import Table
from rich.console import Console
from threading import Lock
from enum import Enum
from dataclasses import dataclass


class TaskStatus(Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


@dataclass
class Task:
    name: str
    total: int = 100
    completed: int = 0
    status: TaskStatus = TaskStatus.PENDING


class TaskDisplay:
    def __init__(self):
        self.tasks: dict[str, Task] = {}
        self.console = Console()
        self._live = None

    def add_task(self, task_id: str, name: str, total: int = 100) -> str:
        self.tasks[task_id] = Task(name=name, total=total)
        return task_id

    def start_task(self, task_id: str):
        self.tasks[task_id].status = TaskStatus.RUNNING
        self._refresh()

    def advance(self, task_id: str, amount: int = 1):
        self.tasks[task_id].completed += amount
        self._refresh()

    def complete(self, task_id: str):
        task = self.tasks[task_id]
        task.completed = task.total
        task.status = TaskStatus.COMPLETED
        self._refresh()

    def fail(self, task_id: str):
        self.tasks[task_id].status = TaskStatus.FAILED
        self._refresh()

    def _refresh(self):
        if self._live:
            self._live.update(self._render())

    def _render(self) -> Table:
        table = Table(show_header=False, box=None, padding=(0, 1))
        table.add_column("Status", width=2)
        table.add_column("Task", style="white", width=70)
        table.add_column("Progress", width=25)

        icons = {
            TaskStatus.PENDING: "[dim]○[/dim]",
            TaskStatus.RUNNING: "[yellow]●[/yellow]",
            TaskStatus.COMPLETED: "[green]✓[/green]",
            TaskStatus.FAILED: "[red]✗[/red]",
        }

        for task in self.tasks.values():
            pct = task.completed / task.total if task.total > 0 else 0
            bar_width = 15
            filled = int(bar_width * pct)

            if task.status == TaskStatus.PENDING:
                bar = f"[dim]{'─' * bar_width}[/dim]"
            elif task.status == TaskStatus.COMPLETED:
                bar = f"[green]{'━' * bar_width}[/green]"
            elif task.status == TaskStatus.FAILED:
                bar = f"[red]{'━' * filled}[/red][dim]{'─' * (bar_width - filled)}[/dim]"
            else:
                bar = f"[yellow]{'━' * filled}[/yellow][dim]{'─' * (bar_width - filled)}[/dim]"

            progress = f"{bar} [dim]{pct:>5.1%}[/dim]"
            name_style = "dim" if task.status == TaskStatus.PENDING else "white"

            table.add_row(icons[task.status], f"[{name_style}]{task.name}[/{name_style}]", progress)

        return table

    def __enter__(self):
        self._live = Live(self._render(), console=self.console, refresh_per_second=10)
        self._live.__enter__()
        return self

    def __exit__(self, *args):
        self._live.__exit__(*args)
        self._live = None
